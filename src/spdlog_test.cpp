#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

int main() 
{
    // Set the default logger to file logger
    auto file_logger = spdlog::basic_logger_mt("basic_logger", "logs/basic.txt");
    spdlog::set_default_logger(file_logger); 
    file_logger->info("Welcome to spdlog!");
    file_logger->flush();
    file_logger->error("Some error message with arg: {}", 1);
    
    file_logger->warn("Easy padding in numbers like {:08d}", 12);
    file_logger->critical("Support for int: {0:d};  hex: {0:x};  oct: {0:o}; bin: {0:b}", 42);
    file_logger->info("Support for floats {:03.2f}", 1.23456);
    file_logger->info("Positional args are {1} {0}..", "too", "supported");
    file_logger->info("{:<30}", "left aligned");
    
    file_logger->debug("This message should be displayed..");    
    
    // change log pattern
    file_logger->set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");
    
    // Compile time log levels
    // define SPDLOG_ACTIVE_LEVEL to desired level
    SPDLOG_TRACE("Some trace message with param {}", 42);
    SPDLOG_DEBUG("Some debug message");
    
               
}