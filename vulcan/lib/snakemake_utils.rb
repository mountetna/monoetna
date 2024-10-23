class Vulcan
  class Snakemake
    module Utils 
      def find_affected_downstream_jobs(dag, scheduled_jobs)
        # This function is supposed to find all jobs that are downstream of the scheduled jobs.
        # This could also be called: elements_after_last_match.
        # It takes in a sorted dag, and a list of jobs. 
        # It returns the elements in the dag that appear after the last job in the list.
        # Recall that the way scheduled_jobs is determined, is by sending all the params and available files
        # to snakemake, and it determines what jobs need to be run. This implies 2 things:
        # 1. jobs downstream from this, need more params or available files to be run.
        # 2. some downstream jobs don't need to re-run because their files or params have not changed. 
        # In either case we just return the downstream jobs. 

        # Recall that the dag is flattened to a list, and it could be the case that not all jobs
        # have a dependency between them. 


        last_matched_index = dag.rindex { |element| scheduled_jobs.include?(element) }
        raise "Cannot find any matching jobs in the dag" if last_matched_index.nil?
  
        # Retrieve elements after the last matched element
        dag[(last_matched_index + 1)..-1]
      end

      module_function :find_affected_downstream_jobs
    end
  end
end
