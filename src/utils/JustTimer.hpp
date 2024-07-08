#pragma once

#include <chrono>
#include <map>
#include <stdexcept>
#include <string_view>
#include <iomanip>
#include <sstream>
#include <iostream>

#ifndef JUST_TIMER_DEFAULT_INSTANCE_ID
    #define JUST_TIMER_DEFAULT_INSTANCE_ID 0
#endif

/********************************************************************
 * @brief Classes for execution time measurement.
 * 
 * Macros for convenient time measurement:
 * 1. Start named timer:
 *   - *JUST_TIMER_START(TimerName)* - Start timer by name. Default instance of GlobalTimer.
 *   - *JUST_TIMER_START_(InstanceID, TimerName)* - Start timer by name for prescribed InstanceID.
 * 
 * 2. Stop named timer:
 *   - *JUST_TIMER_STOP(TimerName)* - Stop timer by name. Default instance of GlobalTimer.
 *   - *JUST_TIMER_STOP_(InstanceID, TimerName)* - Stop timer by name for prescribed InstanceID.
 * 
 * 3. Print timer summary:
 *   - *JUST_TIMER_PRINT(stream)* - Print timer summary to stream. Default instance of GlobalTimer.
 *   - *JUST_TIMER_PRINT_(InstanceID, stream)* - Print timer summary to stream for prescribed InstanceID.
 *********************************************************************/
namespace JustTimer {

    /********************************************************************
     * @brief Wrapper for the set of named timers.
     * @tparam Resolution Timer resolution. Default is std::chrono::milliseconds.
     *********************************************************************/
    template<class Resolution = std::chrono::milliseconds>
    class Timer {
    public:
        using Clock = std::conditional_t<std::chrono::high_resolution_clock::is_steady,
                                        std::chrono::high_resolution_clock,
                                        std::chrono::steady_clock>;

    private:
        // Dummy to avoid gcc bug
        // https://stackoverflow.com/questions/49707184/explicit-specialization-in-non-namespace-scope-does-not-compile-in-gcc
        template<class T, class Dummy = void>
        struct default_suffix { static constexpr std::string_view value() { return ""; } };
        template<class Dummy>
        struct default_suffix<std::chrono::nanoseconds, Dummy> { static constexpr std::string_view value() { return "(ns)"; } };
        template<class Dummy>
        struct default_suffix<std::chrono::microseconds, Dummy> { static constexpr std::string_view value() { return "(µs)"; } };
        template<class Dummy>
        struct default_suffix<std::chrono::milliseconds, Dummy> { static constexpr std::string_view value() { return "(ms)"; } };
        template<class Dummy>
        struct default_suffix<std::chrono::seconds, Dummy> { static constexpr std::string_view value() { return "(s)"; } };
        constexpr static std::string_view suffix_time = default_suffix<Resolution>::value();

        struct SingleTimer {
            Clock::time_point time_start;
            Clock::time_point time_end;
            std::chrono::nanoseconds time_sum;
            size_t count = 0;
            bool is_running = false;
        };

        std::map<std::string, SingleTimer> named_timers;

    public:
        Timer() = default;
        ~Timer() = default;

        Timer(const Timer& et) = delete;
        Timer& operator=(const Timer& et) = delete;

        /** @brief Start timer by name. */
        inline void start(const std::string& timer_name) {
            named_timers[timer_name].time_start = Clock::now();
            named_timers.at(timer_name).is_running = true;
        }

        /** @brief Stop timer by name. */
        inline void stop(const std::string& timer_name) {
            Clock::time_point end_point = Clock::now();

            SingleTimer& current_timer = named_timers.at(timer_name);
            if (current_timer.is_running) {
                current_timer.time_end = end_point;
                current_timer.time_sum += current_timer.time_end - current_timer.time_start;
                current_timer.count++;
                current_timer.is_running = false;
            } else {
                throw std::logic_error("Timer " + timer_name + " is not running");
            }
        }

    public:

        /** @brief Calculate total execution time for the given timer name. */
        auto total_time(const std::string& timer_name) {
            const SingleTimer& timer = named_timers.at(timer_name);
            return std::chrono::duration_cast<Resolution>(timer.time_sum).count();
        }

        /** @brief Calculate number of calls for the given timer name. */
        auto count(const std::string& timer_name) {
            return named_timers.at(timer_name).count;
        }

        /** @brief Clean map of timers */
        void clear(){
            named_timers.clear();
        }

    public:

        /** @brief Print pretty summary table to the given stream (default is std::cout). */
        void print(std::ostream& stream = std::cout) {
            int width_name = 45;
            int width_number = 16;

            std::stringstream ss;
            std::string time_header = "Total Time " + std::string(suffix_time);
            ss << std::left << " " << std::setw(width_name)   << "Timer Name" << " | "
                                << std::setw(width_number) << time_header  << " | "
                                << std::setw(width_number) << "Count"      << "\n";
            stream << std::string(ss.str().size(), '-') << "\n" << ss.str() << std::string(ss.str().size(), '-') << "\n";

            for (const auto& timer : named_timers) {
                auto time_sum = std::chrono::duration_cast<Resolution>(timer.second.time_sum).count();
                stream << std::left << " " << std::setw(width_name)    << timer.first        << " | "
                                        << std::setw(width_number)  << time_sum           << " | "
                                        << std::setw(width_number)  << timer.second.count << "\n";
            }
            stream << std::string(ss.str().size(), '-') << "\n";
        }

        
        /** @brief Print pretty summary table to the given stream (default is std::cout) with ratio related with TotalTimerName */
        void print(std::ostream& stream, const std::string& total_timer_name) {
            int width_name = 45;
            int width_number = 16;

            std::stringstream ss;
            std::string time_header = "Total Time " + std::string(suffix_time);
            ss << std::left << " " << std::setw(width_name)   << "Timer Name" << " | "
                                << std::setw(width_number) << time_header  << " | "
                                << std::setw(width_number) << "Ratio"      << " | "
                                << std::setw(width_number) << "Count"      << "\n";
            stream << std::string(ss.str().size(), '-') << "\n" << ss.str() << std::string(ss.str().size(), '-') << "\n";

            auto total_timer = named_timers[total_timer_name];
            auto total_time = std::chrono::duration_cast<Resolution>(total_timer.time_sum).count();
            for (const auto& timer : named_timers) {
                auto time_sum = std::chrono::duration_cast<Resolution>(timer.second.time_sum).count();
                stream << std::left << " " << std::setw(width_name)    << timer.first        << " | "
                                        << std::setw(width_number)  << time_sum           << " | "
                                        << std::setw(width_number)  << double(time_sum) / total_time << " | "
                                        << std::setw(width_number)  << timer.second.count << "\n";
            }
            stream << std::string(ss.str().size(), '-') << "\n";
        }
    };

    /********************************************************************
     * @brief Global scope wrapper for the set of named timers.
     * @tparam InstanceID ID of current timer instance.
     * @tparam Resolution Timer resolution. Default is std::chrono::milliseconds.
     * 
     * This class designed as singleton with template parameter InstanceID.
     * So it's allow one to create different global timers with defferent
     * instance ID.
     *********************************************************************/
    template<int InstanceID, class Resolution = std::chrono::milliseconds>
    class GlobalTimer : public Timer<Resolution> {
    public:
        static GlobalTimer<InstanceID>& getInstance() {
            static GlobalTimer<InstanceID> instance;
            return instance;
        }

        GlobalTimer(const GlobalTimer<InstanceID>&) = delete;
        GlobalTimer<InstanceID>& operator=(const GlobalTimer<InstanceID>&) = delete;

    protected:
        GlobalTimer() = default;
        ~GlobalTimer() = default;
    };

    /********************************************************************
     * @brief Get instance of global timer with prescribed InstanceID.
     * @tparam InstanceID Prescribed ID of timer instance.
     *********************************************************************/
    template<int InstanceID>
    inline GlobalTimer<InstanceID>& get() {
        return GlobalTimer<InstanceID>::getInstance();
    }

    /********************************************************************
     * @brief Get default instance (InstanceID = 0) of global timer.
     *********************************************************************/
    inline GlobalTimer<JUST_TIMER_DEFAULT_INSTANCE_ID>& get() {
        return GlobalTimer<JUST_TIMER_DEFAULT_INSTANCE_ID>::getInstance();
    }
}

#define JUST_TIMER_START_(InstanceID, TimerName) JustTimer::get<InstanceID>().start(#TimerName)
#define JUST_TIMER_START(TimerName) JUST_TIMER_START_(JUST_TIMER_DEFAULT_INSTANCE_ID, TimerName)

#define JUST_TIMER_STOP_(InstanceID, TimerName) JustTimer::get<InstanceID>().stop(#TimerName)
#define JUST_TIMER_STOP(TimerName) JUST_TIMER_STOP_(JUST_TIMER_DEFAULT_INSTANCE_ID, TimerName)

#define JUST_TIMER_PRINT_(InstanceID, stream) JustTimer::get<InstanceID>().print(stream)
#define JUST_TIMER_PRINT(stream) JUST_TIMER_PRINT_(JUST_TIMER_DEFAULT_INSTANCE_ID, stream)

#define JUST_TIMER_CLEAR_(InstanceID) JustTimer::get<InstanceID>().clear()
#define JUST_TIMER_CLEAR JUST_TIMER_CLEAR_(JUST_TIMER_DEFAULT_INSTANCE_ID)

#define JUST_TIMER_PRINT_RATIO(stream, TotalTimerName) JustTimer::get<JUST_TIMER_DEFAULT_INSTANCE_ID>().print(stream, #TotalTimerName)