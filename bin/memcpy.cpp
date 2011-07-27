#include <cstring>

extern "C" void memcpy_(void *dest, void *src, const int &bytes) {
  if (!bytes) return;
  memcpy(dest, src, bytes);
}
