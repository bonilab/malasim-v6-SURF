// Config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <string>

#include "DrugParameters.h"
#include "EpidemiologicalParameters.h"
#include "GenotypeParameters.h"
#include "ImmuneSystemParameterOverrides.h"
#include "ImmuneSystemParameters.h"
#include "ModelSettings.h"
#include "MosquitoParameters.h"
#include "MovementSettings.h"
#include "ParasiteParameters.h"
#include "PopulationDemographic.h"
#include "PopulationEvents.h"
#include "RaptSettings.h"
#include "SeasonalitySettings.h"
#include "SimulationTimeframe.h"
#include "SpatialSettings/SpatialSettings.h"
#include "StrategyParameters.h"
#include "TherapyParameters.h"
#include "TransmissionSettings.h"

class Config {
public:
  // Disallow copy
  Config(const Config &) = delete;
  Config &operator=(const Config &) = delete;

  // Disallow move
  Config(Config &&) = delete;
  Config &operator=(Config &&) = delete;

  // Constructor and Destructor
  Config() = default;
  virtual ~Config() = default;

  // Load configuration from a YAML file
  bool load(const std::string &filename);

  // Reload configuration (useful for dynamic updates)
  void reload();

  // Validate all cross-field validations
  void validate_all_cross_field_validations();

  // Getters for entire configuration structures
  [[nodiscard]] const ModelSettings &get_model_settings() const { return model_settings_; }
  void set_model_settings(const ModelSettings &settings) { model_settings_ = settings; }
  [[nodiscard]] const SimulationTimeframe &get_simulation_timeframe() const {
    return simulation_timeframe_;
  }
  void set_simulation_timeframe(const SimulationTimeframe &timeframe) {
    simulation_timeframe_ = timeframe;
  }

  [[nodiscard]] const TransmissionSettings &get_transmission_settings() const {
    return transmission_settings_;
  }
  [[nodiscard]] const PopulationDemographic &get_population_demographic() const {
    return population_demographic_;
  }

  void set_population_demographic(const PopulationDemographic &demographic) {
    population_demographic_ = demographic;
  }

  [[nodiscard]] const EpidemiologicalParameters &get_epidemiological_parameters() const {
    return epidemiological_parameters_;
  }
  void set_epidemiological_parameters(const EpidemiologicalParameters &parameters) {
    epidemiological_parameters_ = parameters;
  }

  [[nodiscard]] const ParasiteParameters &get_parasite_parameters() const {
    return parasite_parameters_;
  }
  void set_parasite_parameters(const ParasiteParameters &parameters) {
    parasite_parameters_ = parameters;
  }

  [[nodiscard]] SpatialSettings &get_spatial_settings() {
    /* no const here because Spatial Data class will need to access and modify
     * later */
    return spatial_settings_;
  }
  [[nodiscard]] SeasonalitySettings &get_seasonality_settings() { return seasonality_settings_; }
  [[nodiscard]] MovementSettings &get_movement_settings() { return movement_settings_; }

  [[nodiscard]] const ImmuneSystemParameters &get_immune_system_parameters() const {
    return immune_system_parameters_;
  }
  void set_immune_system_parameters(const ImmuneSystemParameters &parameters) {
    immune_system_parameters_ = parameters;
  }

  [[nodiscard]] const ImmuneSystemParameterOverrides &get_version6_pfpr_incidence_calibrations()
      const {
    return version6_pfpr_incidence_calibrations_;
  }
  [[nodiscard]] bool has_version6_pfpr_incidence_calibrations() const {
    return has_version6_pfpr_incidence_calibrations_;
  }

  [[nodiscard]] GenotypeParameters &get_genotype_parameters() { return genotype_parameters_; }
  [[nodiscard]] const DrugParameters &get_drug_parameters() const { return drug_parameters_; }
  [[nodiscard]] const TherapyParameters &get_therapy_parameters() const {
    return therapy_parameters_;
  }
  [[nodiscard]] StrategyParameters &get_strategy_parameters() { return strategy_parameters_; }
  [[nodiscard]] MosquitoParameters &get_mosquito_parameters() { return mosquito_parameters_; }
  [[nodiscard]] PopulationEvents &get_population_events() { return population_events_; }
  [[nodiscard]] RaptSettings &get_rapt_settings() { return rapt_settings_; }

  // Make relevant getters virtual for mocking
  [[nodiscard]] size_t number_of_locations() const;
  [[nodiscard]] int number_of_age_classes() const;
  static size_t number_of_parasite_types();
  [[nodiscard]] int number_of_tracking_days() const;
  [[nodiscard]] const std::vector<int> &age_structure() const;

  // Keep non-virtual if not mocked directly
  std::vector<Spatial::Location> &location_db();
  // std::vector<IStrategy*> &strategy_db();
  // GenotypeDatabase* genotype_db();

private:
  void reset_load_state();
  void parse_configuration(const YAML::Node &config);
  void parse_version6_pfpr_incidence_calibrations(const YAML::Node &config);
  void select_random_immune_system_parameter_calibration_id();
  void apply_selected_immune_system_parameter_calibration_id();
  void log_genotype_configuration() const;
  static void log_chromosome_configuration(const GenotypeParameters::ChromosomeInfo &chromosome);
  static void log_gene_configuration(const GenotypeParameters::GeneInfo &gene);
  static void log_amino_acid_position_configuration(
      const GenotypeParameters::AminoAcidPosition &amino_acid_position);
  void process_configuration(const YAML::Node &config);
  void parse_population_events(const YAML::Node &config);
  void validate_model_settings() const;
  void validate_simulation_timeframe() const;
  void validate_transmission_settings() const;
  void validate_population_demographic() const;
  void validate_spatial_settings();
  void validate_seasonality_settings();
  void validate_movement_settings() const;
  void validate_parasite_parameters() const;
  void validate_immune_system_parameters() const;
  void validate_genotype_parameters() const;
  void validate_drug_parameters() const;
  [[nodiscard]] int validate_therapy_parameters() const;
  void validate_strategy_parameters() const;
  void validate_strategy_database_ids() const;
  void validate_public_private_strategies() const;
  void validate_epidemiological_parameters(int therapy_max_dosing_days) const;
  void validate_epidemiological_treatment_parameters(int therapy_max_dosing_days) const;
  void validate_epidemiological_transmission_parameters() const;
  void validate_mosquito_parameters() const;
  void validate_population_events() const;

  // Template method for getting a field
  template <typename T>
  [[nodiscard]] const T &get_field(const T &field) const {
    return field;
  }

  // Template method for setting a field
  template <typename T>
  void set_field(T &field, const T &value) {
    field = value;
  }

  // Configuration File Path
  std::string config_file_path_;
  // Model *model_;

  ModelSettings model_settings_;
  TransmissionSettings transmission_settings_;
  PopulationDemographic population_demographic_;
  SimulationTimeframe simulation_timeframe_;
  SpatialSettings spatial_settings_;
  SeasonalitySettings seasonality_settings_;
  MovementSettings movement_settings_;
  ParasiteParameters parasite_parameters_;
  ImmuneSystemParameters immune_system_parameters_;
  ImmuneSystemParameterOverrides version6_pfpr_incidence_calibrations_;
  bool has_version6_pfpr_incidence_calibrations_ = false;
  GenotypeParameters genotype_parameters_;
  DrugParameters drug_parameters_;
  TherapyParameters therapy_parameters_;
  StrategyParameters strategy_parameters_;
  EpidemiologicalParameters epidemiological_parameters_;
  MosquitoParameters mosquito_parameters_;
  PopulationEvents population_events_;
  RaptSettings rapt_settings_;
};

#endif  // CONFIG_H
