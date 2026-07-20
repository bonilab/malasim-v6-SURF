#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>
#include <iomanip>
#include <sstream>
#include "Configuration/ImmuneSystemParameterOverrides.h"

namespace P = ImmuneSystemOverridePaths;

// ---------------------------------------------------------------------------
// Helper: build one calibration_id YAML node by parsing a YAML string.
// This avoids yaml-cpp's colon-path expansion when const char* keys that
// contain colons are used with operator[].  Quoting the keys in the YAML
// text forces yaml-cpp to store them as plain string scalars.
// ---------------------------------------------------------------------------
static YAML::Node make_calibration_id_yaml_node(double p_ci_symp, double z, double kappa,
                                            double midpoint, double p_seek_base,
                                            bool include_mutation_prob = false,
                                            double mutation_prob = -1.0,
                                            bool include_cnv_mult = false,
                                            double cnv_mult = 0.0) {
  std::ostringstream ss;
  ss << std::setprecision(15)
     << "\"" << P::K_P_CI_SYMP   << "\": " << p_ci_symp   << "\n"
     << "\"" << P::K_Z           << "\": " << z            << "\n"
     << "\"" << P::K_KAPPA       << "\": " << kappa        << "\n"
     << "\"" << P::K_MIDPOINT    << "\": " << midpoint     << "\n"
     << "\"" << P::K_P_SEEK_BASE << "\": " << p_seek_base  << "\n";
  if (include_mutation_prob) {
    ss << "\"" << P::K_MUTATION_PROB << "\": " << mutation_prob << "\n";
  }
  if (include_cnv_mult) {
    ss << "\"" << P::K_DEFAULT_CNV_REVERSION_MULTIPLIER << "\": " << cnv_mult << "\n";
  }
  return YAML::Load(ss.str());
}

class YamlImmuneSystemParameterOverridesTest : public ::testing::Test {
protected:
  // Builds a standard 3-calibration_id node (keys 0, 1, 4) using full-path keys.
  // Both genotype override keys are omitted by default (absent from map).
  YAML::Node make_calibration_ids_node(int chosen_calibration_id = 0, bool random_selection = false,
                                  bool include_mutation_prob = false,
                                  double mutation_prob = -1.0,
                                  bool include_cnv_mult = false,
                                  double cnv_mult = 0.0) {
    YAML::Node node;
    node["chosen_calibration_id"] = chosen_calibration_id;
    node["random_selection"]   = random_selection;
    YAML::Node cmap;
    cmap[0] = make_calibration_id_yaml_node(0.5, 3.0, 0.1, 0.2,  0.8,
                                       include_mutation_prob, mutation_prob,
                                       include_cnv_mult,      cnv_mult);
    cmap[1] = make_calibration_id_yaml_node(0.6, 3.0, 0.1, 0.15, 0.6,
                                       include_mutation_prob, mutation_prob,
                                       include_cnv_mult,      cnv_mult);
    cmap[4] = make_calibration_id_yaml_node(1.0, 1.5, 0.5, 0.1,  0.65,
                                       include_mutation_prob, mutation_prob,
                                       include_cnv_mult,      cnv_mult);
    node["calibration_ids"] = cmap;
    return node;
  }
};

// Decode a valid calibration_ids block and check values
TEST_F(YamlImmuneSystemParameterOverridesTest, DecodesValidcalibration_ids) {
  auto node = make_calibration_ids_node(1);
  ImmuneSystemParameterOverrides result;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result));

  EXPECT_EQ(result.get_chosen_calibration_id(), 1);
  EXPECT_FALSE(result.get_random_selection());
  EXPECT_EQ(result.get_calibration_ids().size(), 3u);

  const auto &c0 = result.get_calibration_ids().at(0);
  EXPECT_DOUBLE_EQ(c0.get(P::K_P_CI_SYMP, 0.0),   0.5);
  EXPECT_DOUBLE_EQ(c0.get(P::K_Z, 0.0),           3.0);
  EXPECT_DOUBLE_EQ(c0.get(P::K_KAPPA, 0.0),       0.1);
  EXPECT_DOUBLE_EQ(c0.get(P::K_MIDPOINT, 0.0),    0.2);
  EXPECT_DOUBLE_EQ(c0.get(P::K_P_SEEK_BASE, 0.0), 0.8);
  // genotype keys absent -> not present in map
  EXPECT_FALSE(c0.has(P::K_MUTATION_PROB));
  EXPECT_FALSE(c0.has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER));

  const auto &c1 = result.get_calibration_ids().at(1);
  EXPECT_DOUBLE_EQ(c1.get(P::K_P_CI_SYMP, 0.0),   0.6);
  EXPECT_DOUBLE_EQ(c1.get(P::K_P_SEEK_BASE, 0.0), 0.6);
  EXPECT_FALSE(c1.has(P::K_MUTATION_PROB));
  EXPECT_FALSE(c1.has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER));

  const auto &c4 = result.get_calibration_ids().at(4);
  EXPECT_DOUBLE_EQ(c4.get(P::K_P_CI_SYMP, 0.0),   1.0);
  EXPECT_DOUBLE_EQ(c4.get(P::K_P_SEEK_BASE, 0.0), 0.65);
  EXPECT_FALSE(c4.has(P::K_MUTATION_PROB));
  EXPECT_FALSE(c4.has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER));
}

// Decodes calibration_ids that include both genotype override keys (matches sample_inputs/input.yml)
TEST_F(YamlImmuneSystemParameterOverridesTest, DecodesFullGenotypeOverrides) {
  auto node = make_calibration_ids_node(0, false,
                                   /*include_mutation_prob=*/true, /*mutation_prob=*/0.00085,
                                   /*include_cnv_mult=*/true,      /*cnv_mult=*/0.1);
  ImmuneSystemParameterOverrides result;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result));

  for (const auto &[idx, c] : result.get_calibration_ids()) {
    EXPECT_TRUE(c.has(P::K_MUTATION_PROB))
        << "calibration_id[" << idx << "] should have mutation_probability_per_locus";
    EXPECT_DOUBLE_EQ(c.get(P::K_MUTATION_PROB, 0.0), 0.00085)
        << "calibration_id[" << idx << "] mutation_probability_per_locus";
    EXPECT_TRUE(c.has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER))
        << "calibration_id[" << idx << "] should have default_cnv_reversion_multiplier";
    EXPECT_DOUBLE_EQ(c.get(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER, -1.0), 0.1)
        << "calibration_id[" << idx << "] default_cnv_reversion_multiplier";
  }
}

// has_selected_calibration_id returns true for a present index
TEST_F(YamlImmuneSystemParameterOverridesTest, HasSelectedcalibration_idTrue) {
  auto node = make_calibration_ids_node(0);
  ImmuneSystemParameterOverrides result;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result);

  EXPECT_TRUE(result.has_selected_calibration_id());
  EXPECT_NO_THROW({
    const auto &c = result.get_selected_calibration_id();
    EXPECT_DOUBLE_EQ(c.get(P::K_Z, 0.0), 3.0);
  });
}

// has_selected_calibration_id returns false when chosen_calibration_id is missing
TEST_F(YamlImmuneSystemParameterOverridesTest, HasSelectedcalibration_idFalseWhenMissing) {
  auto node = make_calibration_ids_node(99);  // index 99 not in calibration_ids
  ImmuneSystemParameterOverrides result;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result);

  EXPECT_FALSE(result.has_selected_calibration_id());
  EXPECT_THROW(
      {
        const auto &selected_calibration_id = result.get_selected_calibration_id();
        (void)selected_calibration_id;
      },
      std::runtime_error);
}

// Encode round-trips back to decodable YAML (with both genotype keys)
TEST_F(YamlImmuneSystemParameterOverridesTest, EncodeDecodeRoundtrip) {
  auto node = make_calibration_ids_node(0, false,
                                   /*include_mutation_prob=*/true, /*mutation_prob=*/0.00085,
                                   /*include_cnv_mult=*/true,      /*cnv_mult=*/0.0);
  ImmuneSystemParameterOverrides original;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, original);

  YAML::Node encoded = YAML::convert<ImmuneSystemParameterOverrides>::encode(original);

  ImmuneSystemParameterOverrides decoded;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(encoded, decoded));

  EXPECT_EQ(decoded.get_chosen_calibration_id(), 0);
  EXPECT_FALSE(decoded.get_random_selection());
  EXPECT_EQ(decoded.get_calibration_ids().size(), original.get_calibration_ids().size());
  EXPECT_DOUBLE_EQ(decoded.get_calibration_ids().at(0).get(P::K_P_CI_SYMP, 0.0),
                   original.get_calibration_ids().at(0).get(P::K_P_CI_SYMP, 0.0));
  EXPECT_DOUBLE_EQ(decoded.get_calibration_ids().at(4).get(P::K_P_SEEK_BASE, 0.0),
                   original.get_calibration_ids().at(4).get(P::K_P_SEEK_BASE, 0.0));
  EXPECT_DOUBLE_EQ(decoded.get_calibration_ids().at(0).get(P::K_MUTATION_PROB, 0.0), 0.00085);
  EXPECT_DOUBLE_EQ(decoded.get_calibration_ids().at(0).get(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER, -1.0),
                   0.0);
}

// random_selection is parsed and preserved by encode/decode
TEST_F(YamlImmuneSystemParameterOverridesTest, RandomSelectionRoundtrip) {
  auto node = make_calibration_ids_node(0, true);
  ImmuneSystemParameterOverrides original;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, original);

  EXPECT_TRUE(original.get_random_selection());

  YAML::Node encoded = YAML::convert<ImmuneSystemParameterOverrides>::encode(original);

  ImmuneSystemParameterOverrides decoded;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(encoded, decoded));
  EXPECT_TRUE(decoded.get_random_selection());
}

// random_selection defaults to false when omitted
TEST_F(YamlImmuneSystemParameterOverridesTest, RandomSelectionDefaultsToFalseWhenMissing) {
  auto node = make_calibration_ids_node(0);
  node.remove("random_selection");

  ImmuneSystemParameterOverrides decoded;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(node, decoded));
  EXPECT_FALSE(decoded.get_random_selection());
}

// Missing 'chosen_calibration_id' throws
TEST_F(YamlImmuneSystemParameterOverridesTest, MissingUsedInSimulationThrows) {
  YAML::Node node;
  node["calibration_ids"] = YAML::Node(YAML::NodeType::Map);
  ImmuneSystemParameterOverrides result;
  EXPECT_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result),
               std::runtime_error);
}

// Missing 'calibration_ids' throws
TEST_F(YamlImmuneSystemParameterOverridesTest, Missingcalibration_idsThrows) {
  YAML::Node node;
  node["chosen_calibration_id"] = 0;
  ImmuneSystemParameterOverrides result;
  EXPECT_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result),
               std::runtime_error);
}

// A calibration_id with only some keys is decoded successfully (all keys are optional in the map)
TEST_F(YamlImmuneSystemParameterOverridesTest, calibration_idWithPartialKeysDecodesOk) {
  // Build the whole node from a single YAML document so that `calibration_ids` is a
  // MAP with integer key 0.  Assigning cmap[0] with 0 as the only/first index
  // makes yaml-cpp create a SEQUENCE instead of a map, which then fails to
  // decode (map iterator vs sequence iterator).  Quoting the colon-containing
  // key keeps it stored as a plain scalar string.
  std::ostringstream ss;
  ss << "chosen_calibration_id: 0\n"
     << "calibration_ids:\n"
     << "  0:\n"
     << "    \"" << P::K_P_CI_SYMP << "\": 0.5\n";
  YAML::Node node = YAML::Load(ss.str());

  ImmuneSystemParameterOverrides result;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result));
  EXPECT_TRUE(result.get_calibration_ids().at(0).has(P::K_P_CI_SYMP));
  EXPECT_FALSE(result.get_calibration_ids().at(0).has(P::K_Z));
  EXPECT_FALSE(result.get_calibration_ids().at(0).has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER));
}

// Non-sequential keys are decoded correctly and the key-collection logic used by
// Config::load random_selection picks only valid keys.
TEST_F(YamlImmuneSystemParameterOverridesTest, RandomSelectionNonSequentialKeys) {
  YAML::Node node;
  node["chosen_calibration_id"] = 2;
  node["random_selection"]   = true;
  const std::vector<int> inputKeys = {2, 7, 15};
  YAML::Node cmap;
  for (int key : inputKeys) {
    cmap[key] = make_calibration_id_yaml_node(0.5, 3.0, 0.1, 0.2, 0.8);
  }
  node["calibration_ids"] = cmap;

  ImmuneSystemParameterOverrides result;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result));

  EXPECT_TRUE(result.get_random_selection());
  EXPECT_EQ(result.get_calibration_ids().size(), 3u);

  for (int key : inputKeys) {
    EXPECT_TRUE(result.get_calibration_ids().count(key) > 0)
        << "Expected key " << key << " to exist in calibration_ids";
  }

  std::vector<int> actualKeys;
  actualKeys.reserve(result.get_calibration_ids().size());
  for (const auto &[key, val] : result.get_calibration_ids()) { actualKeys.push_back(key); }

  EXPECT_EQ(actualKeys, inputKeys)
      << "calibration_id keys should be sorted and match the input set";

  for (std::size_t pick = 0; pick < actualKeys.size(); ++pick) {
    const int selectedKey = actualKeys[pick];
    result.set_chosen_calibration_id(selectedKey);
    EXPECT_TRUE(result.has_selected_calibration_id())
        << "pick=" << pick << " -> key=" << selectedKey << " should be a valid calibration_id";
    EXPECT_NO_THROW({
      const auto &selected_calibration_id = result.get_selected_calibration_id();
      (void)selected_calibration_id;
    });
  }

  result.set_chosen_calibration_id(99);
  EXPECT_FALSE(result.has_selected_calibration_id());
}

// Genotype override keys: both absent in YAML means neither key is present in map
TEST_F(YamlImmuneSystemParameterOverridesTest, GenotypeKeysAbsentWhenNotInYaml) {
  auto node = make_calibration_ids_node(0);
  ImmuneSystemParameterOverrides result;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result);

  for (const auto &[idx, c] : result.get_calibration_ids()) {
    EXPECT_FALSE(c.has(P::K_MUTATION_PROB))
        << "calibration_id[" << idx << "] should not have mutation_probability_per_locus key";
    EXPECT_FALSE(c.has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER))
        << "calibration_id[" << idx << "] should not have default_cnv_reversion_multiplier key";
  }
}

// mutation_probability_per_locus: explicit -1 in YAML is present in map with value -1
TEST_F(YamlImmuneSystemParameterOverridesTest, MutationProbExplicitMinusOneIsPreserved) {
  auto node = make_calibration_ids_node(0, false, /*include_mutation_prob=*/true, /*mutation_prob=*/-1.0);
  ImmuneSystemParameterOverrides result;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result);

  for (const auto &[idx, c] : result.get_calibration_ids()) {
    EXPECT_TRUE(c.has(P::K_MUTATION_PROB))
        << "calibration_id[" << idx << "] should have mutation_probability_per_locus key";
    EXPECT_DOUBLE_EQ(c.get(P::K_MUTATION_PROB, 0.0), -1.0)
        << "calibration_id[" << idx << "] explicit -1 should be preserved";
  }
}

// mutation_probability_per_locus: positive value is parsed and round-trips
TEST_F(YamlImmuneSystemParameterOverridesTest, MutationProbPositiveValueRoundtrip) {
  const double inputMutProb = 0.005;
  auto node = make_calibration_ids_node(0, false, /*include_mutation_prob=*/true, inputMutProb);

  ImmuneSystemParameterOverrides original;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, original);

  for (const auto &[idx, c] : original.get_calibration_ids()) {
    EXPECT_TRUE(c.has(P::K_MUTATION_PROB));
    EXPECT_DOUBLE_EQ(c.get(P::K_MUTATION_PROB, 0.0), inputMutProb)
        << "calibration_id[" << idx << "] should have mutation_probability_per_locus=" << inputMutProb;
  }

  YAML::Node encoded = YAML::convert<ImmuneSystemParameterOverrides>::encode(original);
  ImmuneSystemParameterOverrides decoded;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(encoded, decoded));

  for (const auto &[idx, c] : decoded.get_calibration_ids()) {
    EXPECT_DOUBLE_EQ(c.get(P::K_MUTATION_PROB, 0.0), inputMutProb)
        << "round-tripped calibration_id[" << idx << "] should preserve mutation_probability_per_locus";
  }
}

// default_cnv_reversion_multiplier: zero value (as in most calibration_ids in input.yml) is parsed and round-trips
TEST_F(YamlImmuneSystemParameterOverridesTest, CnvReversionMultiplierZeroRoundtrip) {
  const double inputCnvMult = 0.0;
  auto node = make_calibration_ids_node(0, false,
                                   /*include_mutation_prob=*/false, -1.0,
                                   /*include_cnv_mult=*/true, inputCnvMult);

  ImmuneSystemParameterOverrides original;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, original);

  for (const auto &[idx, c] : original.get_calibration_ids()) {
    EXPECT_TRUE(c.has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER));
    EXPECT_DOUBLE_EQ(c.get(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER, -1.0), inputCnvMult)
        << "calibration_id[" << idx << "] default_cnv_reversion_multiplier should be 0.0";
  }

  YAML::Node encoded = YAML::convert<ImmuneSystemParameterOverrides>::encode(original);
  ImmuneSystemParameterOverrides decoded;
  EXPECT_NO_THROW(YAML::convert<ImmuneSystemParameterOverrides>::decode(encoded, decoded));

  for (const auto &[idx, c] : decoded.get_calibration_ids()) {
    EXPECT_DOUBLE_EQ(c.get(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER, -1.0), inputCnvMult)
        << "round-tripped calibration_id[" << idx << "] should preserve default_cnv_reversion_multiplier";
  }
}

// default_cnv_reversion_multiplier: non-zero value (as in calibration_id 0 of input.yml) is parsed correctly
TEST_F(YamlImmuneSystemParameterOverridesTest, CnvReversionMultiplierNonZeroValue) {
  auto node = make_calibration_ids_node(0, false,
                                   /*include_mutation_prob=*/true,  /*mutation_prob=*/0.00085,
                                   /*include_cnv_mult=*/true,       /*cnv_mult=*/0.1);

  ImmuneSystemParameterOverrides result;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result);

  const auto &c0 = result.get_calibration_ids().at(0);
  EXPECT_TRUE(c0.has(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER));
  EXPECT_DOUBLE_EQ(c0.get(P::K_DEFAULT_CNV_REVERSION_MULTIPLIER, -1.0), 0.1);
  EXPECT_TRUE(c0.has(P::K_MUTATION_PROB));
  EXPECT_DOUBLE_EQ(c0.get(P::K_MUTATION_PROB, 0.0), 0.00085);
}

// log_all does not crash (smoke test)
TEST_F(YamlImmuneSystemParameterOverridesTest, LogAllDoesNotCrash) {
  auto node = make_calibration_ids_node(0, false,
                                   /*include_mutation_prob=*/true,  /*mutation_prob=*/0.00085,
                                   /*include_cnv_mult=*/true,       /*cnv_mult=*/0.0);
  ImmuneSystemParameterOverrides result;
  YAML::convert<ImmuneSystemParameterOverrides>::decode(node, result);
  EXPECT_NO_THROW(result.log_all());
}