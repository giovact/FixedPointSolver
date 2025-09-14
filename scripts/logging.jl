using Logging
using Printf  
using Dates

get_current_datetime() =  Dates.format(Dates.now(), "yyyy-u-dd HH:MM:SS") 

# Global lock for ensuring thread-safe file access
const log_lock = ReentrantLock()

create_subdir("logs","")
critical_line_log_file = "logs/critical_line_log.txt"
parallel_scans_log_file = "logs/parallel_scans_log.txt"
exploration_log_file = "logs/exploration_log.txt"

end_log(filepath::String)  = log_message("\n\n\n\n\n",filepath)    


# Function to log messages to a shared file (thread-safe)
function log_message(message::String, file_path::String)
    lock(log_lock)  # Ensure only one thread writes at a time
    open(file_path, "a") do io  # Open the file in append mode
        @printf(io, "%s: %s\n", get_current_datetime(), message)
    end
    unlock(log_lock)
end


#add at the end log_message("All tasks completed", log_file)


function print_and_log(s::String,filepath::String)
    println(s)
    log_message(s,filepath)
end



############### current commit and relative info ##############
using LibGit2

function get_commit_info()
    commit_hash = readchomp(`git rev-parse HEAD`)
    commit_date = readchomp(`git show -s --format=%ci HEAD`)
    commit_author = readchomp(`git show -s --format=%an HEAD`)
    commit_message = readchomp(`git show -s --format=%s HEAD`)
    return Dict("commit_hash"=>commit_hash,
                "commit_date"=>commit_date,
                "commit_author"=>commit_author,
                "commit_message"=>commit_message)

end
