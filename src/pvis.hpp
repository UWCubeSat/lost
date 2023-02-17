#ifndef PVIS_H_
#define PVIS_H_

#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>   // setfill, setw
#include <iostream>  // cout, hex
#include <map>
#include <memory>  // unique_ptr
// #include <stack>
#include <string>
#include <vector>

class Wrapper;

template <class T>
class TypedWrapper;

template <class T>
class TypedWrapperN;

template <class T>
class TypedWrapperVB;

template <class T>
class TypedWrapperCM;

template <class T>
void PrintBytes(T data) {
  int size = sizeof(data);
  std::cout << "Size: " << size << std::endl;
  uint8_t* p = reinterpret_cast<uint8_t*>(&data);

  // Little endian on x86-64
  for (int i = 0; i < size; i++) {
    // must convert byte to unsigned int first
    std::cout << std::hex << std::setfill('0') << std::setw(2) << +*p << " " << std::dec
              << std::flush;
    p += 1;
  }
  std::cout << std::endl;
}

// TODO: version of above function where bytes are added to a vector
// Then reading bytes can print out original String (maybe another tool entirely??)

class StckFrame {
 public:
  std::string funcName;
  std::map<std::string, Wrapper*> localVars;

  StckFrame(std::string funcName = "") : funcName(funcName) {}
};

#define INIT_LOGGER(fname) std::ofstream loggerFOut(fname, std::ios::out);

#define GET_LOGGER extern std::ofstream loggerFOut;

#define CLOSE_LOGGER loggerFOut.close();

// Emulate "collection of different types" by having map to pointer
// typedef std::stack<std::map<std::string, Wrapper*>> fstck;
typedef std::vector<StckFrame> fstck;

static fstck funcStack;

// Neat macro to get the stringified version of name rather than its actual expression
#define VAR_NAME(name) #name

#define TRACK_VAR(varName) \
  TypedWrapperN<decltype(varName)> tw_##varName(&varName, VAR_NAME(varName), funcStack);

#define TRACK_VAR_VB(varName) \
  TypedWrapperVB<decltype(varName)> tw_##varName(&varName, VAR_NAME(varName), funcStack);

#define TRACK_VAR_CM(varName) \
  TypedWrapperCM<decltype(varName)> tw_##varName(&varName, VAR_NAME(varName), funcStack);

class Wrapper {
 public:
  Wrapper(){};

  virtual ~Wrapper() {}

  // "virtual friend" trick
  friend std::ostream& operator<<(std::ostream& out, const Wrapper& wrp) {
    wrp.Print(out);
    wrp.VPrint(out);
    wrp.CMPrint(out);
    return out;
  }

 protected:
  virtual void Print(std::ostream& out) const = 0;
  virtual void VPrint(std::ostream& out) const = 0;
  virtual void CMPrint(std::ostream& out) const = 0;
  // virtual void Write(std::ofstream& f) const = 0;
};

template <class T>
class TypedWrapper : public Wrapper {
 protected:
  T* obj;
  std::string varName;
  fstck& fs;

 public:
  TypedWrapper(T* obj, std::string varName, fstck& fs)
      : Wrapper(), obj(obj), varName(varName), fs(fs) {
    // std::cout << "constructed var: " << varName << std::endl;
    fs.back().localVars[varName] = this;
  }

  // disable default ctor, cctor, and copy-assignment
  TypedWrapper<T>() = delete;
  TypedWrapper<T>(const TypedWrapper& other) = delete;
  TypedWrapper& operator=(const TypedWrapper& other) = delete;

  ~TypedWrapper<T>() {
    // int retval = fs.back().localVars.erase(varName);
    std::cout << "Destroying wrapper: " << varName << std::endl;
    // assert(retval == 1);
  }

  virtual void Print(std::ostream& out) const {}
  virtual void VPrint(std::ostream& out) const {}
  virtual void CMPrint(std::ostream& out) const {}
};

template <class T>
class TypedWrapperN : public TypedWrapper<T> {
 public:
  TypedWrapperN(T* obj, std::string varName, fstck& fs) : TypedWrapper<T>(obj, varName, fs) {}

  virtual void Print(std::ostream& out) const { out << *this->obj << std::flush; }
};

template <class T>
class TypedWrapperVB : public TypedWrapperN<T> {
 public:
  TypedWrapperVB(T* obj, std::string varName, fstck& fs) : TypedWrapperN<T>(obj, varName, fs) {}

  virtual void VPrint(std::ostream& out) const {
    if (*this->obj == nullptr) {
      out << " -> NULL" << std::flush;
    } else
      out << " -> " << *(*this->obj) << std::flush;
  }
};

template <class T>
class TypedWrapperCM : public TypedWrapper<T> {
 public:
  TypedWrapperCM(T* obj, std::string varName, fstck& fs) : TypedWrapper<T>(obj, varName, fs) {}

  virtual void Print(std::ostream& out) const {}

  virtual void CMPrint(std::ostream& out) const {
    int size = sizeof(T);
    out << "(size = " << size << " B)" << std::endl;
    uint8_t* p = reinterpret_cast<uint8_t*>(this->obj);

    // Little endian on x86-64
    for (int i = 0; i < size; i++) {
      // must convert to unsigned int first
      out << std::hex << std::setfill('0') << std::setw(2) << +*p << " " << std::dec << std::flush;
      p += 1;
    }
    out << std::endl;
  }
};

#ifndef NPRINT
#define PRINT_LOCALS(ind)                                                     \
  if (!funcStack.empty()) {                                                   \
    std::cout << "--------------- " << funcStack.back().funcName << "[" << ind \
              << "] ---------------" << std::endl;                            \
    for (const auto& e : funcStack.back().localVars) {                         \
      std::cout << e.first << ": " << *(e.second) << std::endl;               \
    }                                                                         \
    std::cout << "------------------" << std::endl;                           \
  }
#else
#define PRINT_LOCALS(ind)
#endif  // NPRINT

#define LOG(ind)                                                               \
  if (!funcStack.empty()) {                                                    \
    loggerFOut << "--------------- " << funcStack.back().funcName << "[" << ind \
               << "] ---------------" << std::endl;                            \
    for (const auto& e : funcStack.back().localVars) {                          \
      loggerFOut << e.first << ": " << *(e.second) << std::endl;               \
    }                                                                          \
    loggerFOut << "------------------" << std::endl;                           \
  }

//////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_STACK_FRAME_N(funcName) \
  StckFrame localVars(funcName);     \
  funcStack.push_back(localVars);         \

#define MAKE_STACK_FRAME     \
  StckFrame localVars;       \
  funcStack.push_back(localVars); \

#define DEL_STACK_FRAME  funcStack.pop_back();

#endif  // PVIS_H_