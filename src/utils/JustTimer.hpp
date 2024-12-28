#pragma once

#include <chrono>
#include <stdexcept>
#include <string_view>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fmt/core.h>
#include <fmt/chrono.h>
#include <fmt/args.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/compile.h>

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

        struct SingleTimer{
            std::string name;
            Clock::time_point time_start;
            Clock::time_point time_end;
            std::chrono::nanoseconds time_sum = std::chrono::nanoseconds(0);
            size_t count = 0;
            bool is_running = false;
            int level = 0; // For nested timers (to beatify printing)
            SingleTimer(const std::string& name_, int level_) : name(name_), level(level_){}
        };

        std::vector<SingleTimer> timers;
        int level = -1;
        int max_level = -1;

    public:
        Timer() = default;
        ~Timer() = default;

        Timer(const Timer& et) = delete;
        Timer& operator=(const Timer& et) = delete;

        inline bool is_new_timer(const std::string& timer_name){
            if (timers.size() == 0) return true;
            auto it = std::find_if(timers.begin(), timers.end(),
                [&timer_name](const SingleTimer& x){
                    return x.name == timer_name;
                }
            );
            return it == timers.end();
        }

        SingleTimer& get_timer(const std::string& timer_name){
            auto it = std::find_if(timers.begin(), timers.end(),
                [&timer_name](const SingleTimer& x){
                    return x.name == timer_name;
                }
            );
            return *it;
        }

        /** @brief Start timer by name. */
        inline void start(const std::string& timer_name){
            if (is_new_timer(timer_name)){
                level++;
                timers.push_back(SingleTimer(timer_name, level));
                max_level = std::max(max_level, level);
            }
            SingleTimer& tptr = get_timer(timer_name);
            if (tptr.is_running){
                throw std::logic_error("Timer " + timer_name + " is running! Must be stop before restarting\n");
            }
            tptr.time_start = Clock::now();
            tptr.is_running = true;
        }

        /** @brief Stop timer by name. */
        inline void stop(const std::string& timer_name) {
            Clock::time_point end_point = Clock::now();
            if (is_new_timer(timer_name)){
                throw std::logic_error("Timer " + timer_name + " doesn't exist\n");
            }
            SingleTimer& tptr = get_timer(timer_name);
            if (tptr.is_running){
                tptr.time_end = end_point;
                tptr.time_sum += tptr.time_end - tptr.time_start;
                if (tptr.count == 0) level--;
                tptr.count++;
                tptr.is_running = false;
            } else {
                throw std::logic_error("Timer " + timer_name + " is not running\n");
            }
        }

    public:

        /** @brief Calculate total execution time for the given timer name. */
        auto total_time(const std::string& timer_name) {
            SingleTimer& timer = get_timer(timer_name);
            return std::chrono::duration_cast<Resolution>(timer.time_sum).count();
        }

        /** @brief Calculate number of calls for the given timer name. */
        auto count(const std::string& timer_name) {
            return get_timer(timer_name).count;
        }

        /** @brief Clean map of timers */
        void clear(){
            timers.clear();
            level = -1;
            max_level = -1;
        }

    public:       
        /** @brief Print pretty summary table to the given stream (default is std::cout) with ratio related with TotalTimerName */
        void print(std::ostream& stream){
            int level_tab = 2;
            int level_tabs = std::max(0, level_tab * max_level);
            int width_name = 30 + level_tabs;
            int width_time = 10;
            int width_ratio = 7;
            int width_count = 7;
            std::string time_header = "Time " + std::string(suffix_time);
            std::string line = fmt::format(" {:<{}} | {:<{}} | {:<{}} | {:<{}} \n", "Timer Name", width_name, time_header, width_time, "Ratio", width_ratio, "Count", width_count);
            stream << std::string(line.size(), '-') << "\n" << line << std::string(line.size(), '-') << "\n";

            const SingleTimer& total_timer = timers.front();
            double total_time = std::chrono::duration_cast<Resolution>(total_timer.time_sum).count();
            for (const auto& timer : timers) {
                double time_sum = std::chrono::duration_cast<Resolution>(timer.time_sum).count();
                std::string timer_name = std::string(timer.level * level_tab, '-') + " " + timer.name;
                std::string lline = fmt::format(" {:<{}} | {:<{}} | {:<{}.3f} | {:<{}} \n", timer_name, width_name, time_sum, width_time, double(time_sum) / total_time, width_ratio, timer.count, width_count);
                stream << lline;
            }
            stream << std::string(line.size(), '-') << "\n";
        }
        // void print(std::ostream& stream){
        //     int level_tab = 2;
        //     int level_tabs = std::max(0, level_tab * max_level);
        //     int width_name = 30 + level_tabs;
        //     constexpr int width_time = 10;
        //     constexpr int width_ratio = 7;
        //     constexpr int width_count = 7;
        //     std::string time_header = "Time " + std::string(suffix_time);
        //     std::stringstream ss;
        //     std::string line = fmt::format(" {:<{}} | {:<{}} | {:<{}} | {:<{}} ", "Timer Name", width_name, time_header, width_time, "Ratio", width_ratio, "Count", width_count);
        //     ss  << std::setw(width_name) << std::left << " " << "Timer Name" << " | "
        //         << std::setw(width_time) << time_header << " | "
        //         << std::setw(width_ratio) << "Ratio" << " | "
        //         << std::setw(width_count) << "Count" << "\n";
        //     stream << std::string(ss.str().size(), '-') << "\n" << ss.str() << std::string(ss.str().size(), '-') << "\n";

        //     const SingleTimer& total_timer = timers.front();
        //     auto total_time = std::chrono::duration_cast<Resolution>(total_timer.time_sum).count();
        //     for (const auto& timer : timers) {
        //         auto time_sum = std::chrono::duration_cast<Resolution>(timer.time_sum).count();
        //         std::string timer_name = std::string(timer.level * level_tab, '-') + " " + timer.name;
        //         stream  << std::left << " " << std::setw(width_name) << timer_name << " | "
        //                 << std::setw(width_time) << std::setprecision(std::min(15, width_time)) << time_sum << " | "
        //                 << std::setw(width_ratio) << std::setprecision(std::min(15, width_ratio)) << double(time_sum) / total_time << " | "
        //                 << std::setw(width_count) << std::setprecision(std::min(15, width_count)) << timer.count << "\n";
        //     }
        //     stream << std::string(ss.str().size(), '-') << "\n";
        // }
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