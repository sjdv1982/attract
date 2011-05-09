#include <cstring>

extern "C" void memcpy_(void *dest, void *src, const int &bytes) {
  memcpy(dest, src, bytes);
}
