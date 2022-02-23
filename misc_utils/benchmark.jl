function benchmark()
   println("======================= First run:")
   @time main()

   println("\n\n======================= Second run:")
   Profile.init(delay=0.1)
   Profile.clear()
   Profile.clear_malloc_data()
   @profile @time main()

   r = Profile.retrieve()
   f = open("output_file_name.txt", "w")
   #Profile.print(f; combine=true, mincount=100)
   Profile.print(IOContext(f, :displaysize => (200,500)); mincount=100,combine=true)
   close(f)
end
