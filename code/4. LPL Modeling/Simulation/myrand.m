function [ output ] = myrand(  )



global nearest_time nearest_time_last rand_save;
global DEBUG_RNG;


    if DEBUG_RNG
        if nearest_time == nearest_time_last
            output = rand_save;
        else
            nearest_time_last = nearest_time;
            rand_save = rand;
            output = rand_save;
        end
    else
        output = rand;
    end

end

