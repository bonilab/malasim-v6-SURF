#ifndef ICONFIG_DATA_H
#define ICONFIG_DATA_H

class IConfigData {
public:
  IConfigData() = default;

  IConfigData(const IConfigData &) = default;
  IConfigData &operator=(const IConfigData &) = default;

  IConfigData(IConfigData &&) noexcept = default;
  IConfigData &operator=(IConfigData &&) noexcept = default;

  virtual ~IConfigData() = default;
  virtual void process_config() = 0;
};
#endif  // ICONFIG_DATA_H
