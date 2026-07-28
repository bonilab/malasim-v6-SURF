#ifndef ICONFIG_DATA_H
#define ICONFIG_DATA_H

class IConfigData {
public:
  IConfigData(const IConfigData &) = delete;
  IConfigData &operator=(const IConfigData &) = delete;
  IConfigData(IConfigData &&) = delete;
  IConfigData &operator=(IConfigData &&) = delete;

  virtual ~IConfigData() = default;

  virtual void process_config() = 0;

protected:
  IConfigData() = default;
};
#endif  // ICONFIG_DATA_H
