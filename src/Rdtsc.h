#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__i386__) && !defined(__x86_64__) && !defined(__sparc__)
#warning No supported architecture found -- timers will return junk.
#endif

static __inline__ uint64_t curtick() {
	uint64_t tick;
#if defined(__i386__)
	unsigned long lo, hi;
	__asm__ __volatile__ (".byte 0x0f, 0x31" : "=a" (lo), "=d" (hi));
	tick = (uint64_t) hi << 32 | lo;
#elif defined(__x86_64__)
	unsigned long lo, hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	tick = (uint64_t) hi << 32 | lo;
#elif defined(__sparc__)
	__asm__ __volatile__ ("rd %%tick, %0" : "=r" (tick));
#endif
	return tick;
}

static __inline__ void startTimer(uint64_t* t) {
	*t = curtick();
}

static __inline__ void stopTimer(uint64_t* t) {
	*t = curtick() - *t;
}

#ifdef __cplusplus
}
#endif

