using ProgressMeter

#=
# Create a progress meter for an unknown number of steps, using a spinner
p = ProgressUnknown(spinner=true, desc="Processing:")

# Simulate a long-running task
for i in 1:2
    # Do some work
    sleep(0.5)
    # Call next! to advance the spinner and update the display
    next!(p)
end

# Call finish! when the task is complete
finish!(p)


n = 10
p = ProgressUnknown(spinner=true, "Computing:")
x = 1

for iter in 1:n
    # Perform some computation
    global x += 2
    sleep(0.5) 

    # Update the progress bar and display additional information
    next!(p; showvalues = [(:iteration_count, iter), (:current_x_value, x)])
end

finish!(p)
=#

p = ProgressUnknown(desc="Working...", spinner=true)

println("Starting work...")
sleep(4)  # Simulate some work
update!(p)
sleep(4)  # Simulate some more work
update!(p)

finish!(p)