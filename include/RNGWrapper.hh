#ifndef RNGWrapper_h
#define RNGWrapper_h 1

#include "tls.hh" // Chứa macro G4ThreadLocal của Geant4

template <class T> class RNGWrapper {
public:
  static void set(T *object, double (T::*func)(void));
  static double rng(void);

private:
  // Thêm G4ThreadLocal để mỗi Worker Thread có một bản sao độc lập
  static G4ThreadLocal T *m_obj;
  static G4ThreadLocal double (T::*m_func)(void);
};

// Khởi tạo các biến tĩnh cục bộ luồng
template <class T> G4ThreadLocal T *RNGWrapper<T>::m_obj = nullptr;
template <class T> G4ThreadLocal double (T::*RNGWrapper<T>::m_func)(void) = nullptr;

template <class T> void RNGWrapper<T>::set(T *object, double (T::*func)(void)) {
  m_obj = object;
  m_func = func;
}

template <class T> double RNGWrapper<T>::rng(void) {
  return (m_obj->*m_func)();
}

#endif