#include "debug.h"

#include "platform.h"
#include "../core/logging.h"

#include <execinfo.h>
#include <dlfcn.h>
#include <cxxabi.h>
#include <sstream>

namespace carl {

void printStacktrace(bool interaction) {
	std::stringstream cmd;
	cmd << "gdb --pid=" << getpid() << " -ex bt";
	if (!interaction) cmd << " --batch --quiet";
	system(cmd.str().c_str());
}

inline std::string demangle(const char* name) {
	int status = -4;
	std::unique_ptr<char, void(*)(void*)> res {
		abi::__cxa_demangle(name, NULL, NULL, &status),
		std::free
	};
	return (status == 0) ? res.get() : name ;
}

std::string callingFunction() {
	void* frames[3];
    int cnt = backtrace(frames, sizeof(frames) / sizeof(void*));
	if (cnt == 0) return "<unknown, maybe corrupt>";
	char** symbols = backtrace_symbols(frames, cnt);
	
	std::stringstream out;
	Dl_info info;
	if (dladdr(frames[2], &info) && info.dli_sname) {
		out << demangle(info.dli_sname) << std::endl;
	} else {
		out << "??? " << demangle(symbols[2]) << std::endl;
	}
	free(symbols);
	return out.str();
}

std::string last_assertion_string = "";
int last_assertion_code = 23;

#ifndef NDEBUG
/**
 * Actual signal handler.
 */
void handle_signal(int signal) {
	//printStacktrace(false);
	std::cerr << std::endl << "Catched SIGABRT " << signal << ", exiting with " << (last_assertion_code%256) << std::endl;
	if (last_assertion_string.size() != 0) {
		std::cerr << "Last Assertion catched is: " << last_assertion_string << std::endl;
		std::cerr << "Please check if this is the assertion that is actually thrown." << std::endl;
	}
	exit(last_assertion_code % 256);
}
/**
 * Installs the signal handler.
 */
bool install_signal_handler() {
	CARL_LOG_INFO("carl.util", "Installing signal handler for SIGABRT");
	std::signal(SIGABRT, handle_signal);
	return true;
}
/**
 * Static variable that ensures that install_signal_handler is called.
 */
static bool signal_installed = install_signal_handler();
#endif

}
