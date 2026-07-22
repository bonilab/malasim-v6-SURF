#include "Model.h"

#include <Core/Scheduler/Scheduler.h>
#include <Population/Population.h>
#include <Utils/Random.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "Configuration/Config.h"
#include "MDC/ModelDataCollector.h"
#include "Mosquito/Mosquito.h"
#include "Reporters/Reporter.h"
#include "Treatment/LinearTCM.h"
#include "Treatment/SteadyTCM.h"

bool Model::initialize() {
  config_ = std::make_unique<Config>();
  random_ = std::make_unique<utils::Random>(nullptr, -1);
  scheduler_ = std::make_unique<Scheduler>();
  population_ = std::make_unique<Population>();
  mdc_ = std::make_unique<ModelDataCollector>();
  mosquito_ = std::make_unique<Mosquito>();

  progress_to_clinical_update_function_ = std::make_unique<ClinicalUpdateFunction>(this);
  immunity_clearance_update_function_ = std::make_unique<ImmunityClearanceUpdateFunction>(this);
  having_drug_update_function_ = std::make_unique<ImmunityClearanceUpdateFunction>(this);
  clinical_update_function_ = std::make_unique<ImmunityClearanceUpdateFunction>(this);
  reporters_.clear();

  genotype_db_ = std::make_unique<GenotypeDatabase>();
  drug_db_ = std::make_unique<DrugDatabase>();

  if (cli_input_.input_path.empty()) {
    // spdlog::error("Input path is empty. Please provide a valid input path.");
    // return false;
    spdlog::warn("Input path is empty. Model initialized without configuration file.");
    return true;
  }

  // if input path is not empty, load configuration file
  spdlog::info("Loading configuration file: " + cli_input_.input_path);
  if (config_->load(cli_input_.input_path)) {
    if (config_->get_model_settings().get_initial_seed_number() <= 0) {
      random_->set_seed(std::chrono::system_clock::now().time_since_epoch().count());
    } else {
      random_->set_seed(config_->get_model_settings().get_initial_seed_number());
    }

    spdlog::info("Model initialized with seed: " + std::to_string(random_->get_seed()));

    if (cli_input_.output_path.empty()) { cli_input_.output_path = "./"; }

    // add reporter here
    if (cli_input_.reporter.empty()) {
      add_reporter(Reporter::make_report(Reporter::ReportType::SQLITE_MONTHLY_REPORTER));
    } else {
      if (Reporter::report_type_map.contains(cli_input_.reporter)) {
        add_reporter(Reporter::make_report(Reporter::report_type_map[cli_input_.reporter]));
      }
    }

#ifdef ENABLE_TRAVEL_TRACKING
    add_reporter(Reporter::MakeReport(Reporter::TRAVEL_TRACKING_REPORTER));
#endif

    // initialize reporters
    for (auto &reporter : reporters_) {
      reporter->initialize(cli_input_.job_number, cli_input_.output_path);
    }
    spdlog::info("Model initialized reporters.");

    scheduler_->initialize(config_->get_simulation_timeframe().get_starting_date(),
                           config_->get_simulation_timeframe().get_ending_date());
    spdlog::info("Model initialized scheduler.");

    set_treatment_strategy(config_->get_strategy_parameters().get_initial_strategy_id());
    spdlog::info("Model initialized treatment strategy.");

    set_second_line_strategy(config_->get_strategy_parameters().get_second_line_strategy_id());
    spdlog::info("Model initialized second-line treatment strategy.");

    build_initial_treatment_coverage();
    spdlog::info("Model initialized treatment coverage model.");

    mdc_->initialize();
    spdlog::info("Model initialized data collector.");

    spdlog::info("Model initializing population...");
    population_->initialize();
    spdlog::info("Model initialized population.");

    config_->get_movement_settings().get_spatial_model()->prepare();
    spdlog::info("Model initialized movement model.");

    mosquito_->initialize(config_.get());
    spdlog::info("Model initialized mosquito.");

    population_->introduce_initial_cases();
    spdlog::info("Model initialized initial cases.");

    // Take ownership of the events from the config
    auto population_events = config_->get_population_events().release_events();
    for (auto &event : population_events) {
      if (event) {  // Check if the pointer is valid before using
        spdlog::info("Scheduling population event: {} at {}", event->name(), event->get_time());
        scheduler_->schedule_population_event(std::move(event));
      } else {
        spdlog::warn("Encountered a null event pointer during initialization.");
      }
    }

    if (cli_input_.record_movement) {
      // Generate a movement reporter
      auto reporter = Reporter::make_report(Reporter::ReportType::MOVEMENT_REPORTER);
      reporter->initialize(cli_input_.job_number, cli_input_.output_path);
      add_reporter(std::move(reporter));
    }
    is_initialized_ = true;
  } else {
    spdlog::error("Failed to load configuration file: " + cli_input_.input_path);
  }
  return is_initialized_;
}

void Model::release() {
  // Clean up the memory used by the model

  treatment_strategy_ = nullptr;
  second_line_strategy_ = nullptr;
  treatment_coverage_.reset();

  progress_to_clinical_update_function_.reset();
  immunity_clearance_update_function_.reset();
  having_drug_update_function_.reset();
  clinical_update_function_.reset();

  drug_db_.reset();
  genotype_db_.reset();
  mosquito_.reset();
  mdc_.reset();
  population_.reset();
  random_.reset();
  scheduler_.reset();
  config_.reset();

  // simply clear the vector, the unique_ptr will be deleted automatically
  reporters_.clear();
}

void Model::run() {
  if (!is_initialized_) {
    throw std::runtime_error("Model is not initialized. Call Initialize() first.");
  }
  before_run();
  scheduler_->run();
  after_run();
}

void Model::before_run() {
  spdlog::info("Perform before run events");

  // --------------------------------------------------------------------------
  // Verification dump: print every value that was overridden through
  // version6_pfpr_incidence_calibrations for the selected calibration_id, and compare
  // each one against the value that is actually live in the config right now,
  // so we can confirm all overrides were applied correctly before the run.
  // --------------------------------------------------------------------------
  if (config_ != nullptr) {
    namespace isop = ImmuneSystemOverridePaths;

    // NOTE: adjust these two accessor names if your Config exposes the object
    // under a different name (e.g. get_version6_pfpr_incidence_calibrations() /
    // has_immune_system_parameter_calibration_ids()).
    const bool section_present = config_->has_version6_pfpr_incidence_calibrations();
    const auto &overrides = config_->get_version6_pfpr_incidence_calibrations();

    spdlog::info("===== version6_pfpr_incidence_calibrations (before_run verification) =====");
    spdlog::info("  section_present    = {}", section_present);
    spdlog::info("  random_selection   = {}", overrides.get_random_selection());
    spdlog::info("  chosen_calibration_id = {}", overrides.get_chosen_calibration_id());

    if (!section_present) {
      spdlog::info("  No overrides section -> running with default immune-system parameters.");
    } else if (!overrides.has_selected_calibration_id()) {
      spdlog::warn(
          "  chosen_calibration_id={} not found among calibration_ids -> NO overrides applied!",
          overrides.get_chosen_calibration_id());
    } else {
      // Resolve the value currently live in the config for a given override path.
      // Returns NaN for a path with no corresponding live field.
      auto effective_value = [&](const std::string &path) -> double {
        if (path == isop::K_Z) {
          return config_->get_immune_system_parameters()
              .get_immune_effect_on_progression_to_clinical();
        }
        if (path == isop::K_KAPPA) {
          return config_->get_immune_system_parameters().get_factor_effect_age_mature_immunity();
        }
        if (path == isop::K_MIDPOINT) {
          return config_->get_immune_system_parameters().get_midpoint();
        }
        if (path == isop::K_P_CI_SYMP) {
          return config_->get_epidemiological_parameters()
              .get_allow_new_coinfection_to_cause_symptoms()
              .get_probability();
        }
        if (path == isop::K_P_SEEK_BASE) {
          return config_->get_epidemiological_parameters()
              .get_age_based_probability_of_seeking_treatment()
              .get_power()
              .base;
        }
        if (path == isop::K_MUTATION_PROB) {
          return config_->get_genotype_parameters().get_mutation_probability_per_locus();
        }
        if (path == isop::K_DEFAULT_CNV_REVERSION_MULTIPLIER) {
          return config_->get_genotype_parameters().get_default_cnv_reversion_multiplier();
        }
        return std::numeric_limits<double>::quiet_NaN();
      };

      const auto &calibration_id = overrides.get_selected_calibration_id();
      spdlog::info("  selected calibration_id[{}] has {} override(s):",
                   overrides.get_chosen_calibration_id(), calibration_id.overrides.size());

      for (const auto &[path, override_val] : calibration_id.overrides) {
        const double effective = effective_value(path);

        if (std::isnan(effective)) {
          spdlog::warn("    [UNKNOWN PATH] {} = {} (no live config field to verify)", path,
                       override_val);
          continue;
        }

        // mutation_probability_per_locus and default_cnv_reversion_multiplier use
        // a <0 sentinel meaning "keep the existing default" (see
        // Config::apply_selected_immune_system_parameter_calibration_id), so a negative
        // override is intentionally NOT applied.
        const bool sentinel_keep_default =
            ((path == isop::K_MUTATION_PROB) || (path == isop::K_DEFAULT_CNV_REVERSION_MULTIPLIER))
            && (override_val < 0.0);

        if (sentinel_keep_default) {
          spdlog::info("    [KEPT DEFAULT] {}: override={} (<0 sentinel), effective={}", path,
                       override_val, effective);
          continue;
        }

        const bool applied =
            std::fabs(effective - override_val) <= (1e-9 * std::max(1.0, std::fabs(override_val)));
        spdlog::info("    [{}] {}: override={}, effective={}", applied ? "OK" : "MISMATCH", path,
                     override_val, effective);
        if (!applied) { spdlog::warn("    ^ override for '{}' was NOT applied correctly!", path); }
      }
    }
    spdlog::info("======================================================================");
  }

  for (auto &reporter : reporters_) { reporter->before_run(); }
}

void Model::after_run() {
  spdlog::info("Perform after run events");

  mdc_->update_after_run();

  for (auto &reporter : reporters_) { reporter->after_run(); }
}

void Model::begin_time_step() {
  // reset daily variables
  mdc_->begin_time_step();
  report_begin_of_time_step();
}

void Model::end_time_step() {
  // update / calculate daily UTL
  mdc_->end_of_time_step();

  // check to switch strategy
  treatment_strategy_->update_end_of_time_step();
  if (second_line_strategy_ != nullptr && second_line_strategy_ != treatment_strategy_) {
    second_line_strategy_->update_end_of_time_step();
  }

  report_after_time_step();
}

void Model::daily_update() {
  const auto tracking_index = scheduler_->current_time() % config_->number_of_tracking_days();
  const auto circulation_context = Population::prepare_circulation_context();

  for (int location = 0; location < config_->number_of_locations(); ++location) {
    population_->prepare_daily_state_at_location(location);
    population_->perform_infection_event_at_location(location, tracking_index);
    population_->perform_circulation_from_location(location, circulation_context);
    mosquito_->infect_new_cohort_at_location(config_.get(), random_.get(), population_.get(),
                                             location, tracking_index);
    population_->persist_force_of_infection_at_location(location, tracking_index);
  }
}

void Model::monthly_update() {
  monthly_report();

  // reset monthly variables
  mdc_->monthly_update();

  //
  treatment_strategy_->monthly_update();
  if (second_line_strategy_ != nullptr && second_line_strategy_ != treatment_strategy_) {
    second_line_strategy_->monthly_update();
  }

  // update treatment coverage
  treatment_coverage_->monthly_update();
}

void Model::yearly_update() { mdc_->yearly_update(); }

void Model::monthly_report() {
  mdc_->perform_population_statistic();

  for (auto &reporter : reporters_) { reporter->monthly_report(); }
}

void Model::report_begin_of_time_step() {
  for (auto &reporter : reporters_) { reporter->begin_time_step(); }
}

void Model::report_after_time_step() {
  for (auto &reporter : reporters_) { reporter->after_time_step(); }
}

void Model::add_reporter(std::unique_ptr<Reporter> reporter) {
  reporter->set_model(this);
  reporters_.push_back(std::move(reporter));
}

IStrategy* Model::get_treatment_strategy() { return get_instance()->treatment_strategy_; }

void Model::set_treatment_strategy(const int &strategy_id) {
  treatment_strategy_ = strategy_id == -1 ? nullptr : Model::get_strategy_db()[strategy_id].get();
  assert(treatment_strategy_ != nullptr);
  treatment_strategy_->adjust_started_time_point(Model::get_scheduler()->current_time());
}

IStrategy* Model::get_second_line_strategy() { return get_instance()->second_line_strategy_; }

void Model::set_second_line_strategy(const int strategy_id) {
  second_line_strategy_ = strategy_id == -1 ? nullptr : Model::get_strategy_db()[strategy_id].get();
  if (second_line_strategy_ != nullptr && second_line_strategy_ != treatment_strategy_) {
    second_line_strategy_->adjust_started_time_point(Model::get_scheduler()->current_time());
  }
}

ITreatmentCoverageModel* Model::get_treatment_coverage() {
  return get_instance()->treatment_coverage_.get();
}

void Model::set_treatment_coverage(std::unique_ptr<ITreatmentCoverageModel> tcm) {
  if (treatment_coverage_.get() != tcm.get()) {
    if (tcm->p_treatment_under_5.empty() || tcm->p_treatment_over_5.empty()) {
      // copy current value
      tcm->p_treatment_under_5 = treatment_coverage_->p_treatment_under_5;
      tcm->p_treatment_over_5 = treatment_coverage_->p_treatment_over_5;
    }

    if (auto* linear_tcm = dynamic_cast<LinearTCM*>(tcm.get())) {
      linear_tcm->update_rate_of_change();
    }
  }
  treatment_coverage_ = std::move(tcm);
}

void Model::build_initial_treatment_coverage() {
  auto tcm_ptr = std::make_unique<SteadyTCM>();
  for (const auto &location : config_->location_db()) {
    tcm_ptr->p_treatment_under_5.push_back(location.p_treatment_under_5);
    tcm_ptr->p_treatment_over_5.push_back(location.p_treatment_over_5);
  }
  set_treatment_coverage(std::move(tcm_ptr));
}

std::vector<std::unique_ptr<Reporter>> &Model::get_reporters() { return reporters_; }
