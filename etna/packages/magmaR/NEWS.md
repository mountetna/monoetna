# magmaR 1.0.4

* Added 'updateValues()' and 'updateFromDF()' compatibility for table-type models.
* BugFix: Fixed how NULL values are converted to JSON (now as null instead of {}), which allows for data values to be blanked via 'updateValues()' and 'updateFromDF()'.  
* BugFix: Fixed code and documentation for 'retrieveProjects()' after its target endpoint moved to <janus>/api/user/projects.

# magmaR 1.0.3

* Added 'showDisconnected' input to retrieve functions.
* Added 'dryRun' and 'autolink' inputs to update functions.

# magmaR 1.0.2

* Added new data upload function 'updateFromDF()'
* Minor changes to updates-summarization by 'update*()' functions.
* Updated the "Upload" vignette to include this new function.
* Added a `NEWS.md` file to track changes to the package.

# magmaR 1.0.1

* Submitted to CRAN
