#ifndef ICONFIG_DATA_H
#define ICONFIG_DATA_H

class IConfigData {
public:
  IConfigData() = default;

  virtual ~IConfigData() = default;
  virtual void process_config() = 0;
};
#endif  // ICONFIG_DATA_H
