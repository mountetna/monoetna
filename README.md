# Polyphemus

Polyphemus was one of three Cyclops, sons of Poseidon who worked the forges at Mount Etna along with Hephaestus (Vulcan).

Here Polyphemus is a data worker. He browses through [magma](https://github.com/mountetna/magma), which is constantly being updated from various sources, and looks for data which is ready to be run through one of many data loaders. Ideally he does this by consuming a journal of some sort (or being directly notified by other consumers of magma).

The data loaders are attached directly to magma models - they each have their own logic for determining whether they should run. Polyphemus merely runs the checks specified by each loader, which should verify that its inputs exist and are well formed, the sky is turned on, etc. If they are found to pass, and the expected output data is missing (or perhaps out-of-date), Polyphemus will suggest running the loader (or some of several matching loaders).

Polyphemus talks through email, a web interface, or (most likely) hooks into Timur, where he can offer help and use of tools.

A common task Polyphemus might run would be orphan data collection. Nodes found in a Magma project which inappropriately do not have parents (which Polyphemus may discover via occasional polling queries, or notifications from Timur, Magma, etc.) may be suggested for 'orphan data review', during which Polyphemus can delete the orphaned nodes, after making a copy of their data.

However, Polyphemus can coordinate the performance of any task. For example, after uploading some fastq files Polyphemus might suggest running an exome analysis on the fresh input dataset, and can help the user to configure an exome analysis pipeline.

## Managing analysis

Polyphemus publishes, at a web-accessible location, a JSON document describing the analysis. For example:

    {
       analysis: "exome_pipe",
       requester: "dr.who@ucsf.edu",
       name: "dt_lung_tumor_pair",
       request-date: "2016-05-04@12:22",
       status: "new",
       config: {
        project: "dt_lung",
        inputs: [
          "blood_normal.exome",
          "tumor.exome"
       }
    }

The format could be anything, it must only have a unique name, requester, analysis, request date and status.

Valid consumers of Polyphemus (authenticated via HMAC) may view these open requests for analysis and identify available work (however they choose). A simple way to do this might be a cron job that polls every 5 minutes to check Polyphemus's JSON endpoint. If they find something, they respond with a bid for work:

    {
      analysis: "exome_pipe",
      analyst: "saurabh.asthana@ucsf.edu",
      name: "dt_lung_tumor_pair",
      time: "3600"
    }
Polyphemus replies with an acknowledgement to run:

    {
      message: "go"
    }
and updates the analysis description:

    {
      status: "running",
      analyst: "saurabh.asthana@ucsf.edu",
      start-date: "2016-05-04@12:23",
      due-date: "2016-05-04@13:23"
    }

At the expected due date, Polyphemus will run the data loader check again. If it does not pass, Polyphemus will wait 10% longer (extend the due-date). If it still does not pass, Polyphemus will (perhaps, after an even longer wait) defer to the user who requested the analysis, who can either cancel or allow Polyphemus to coordinate help from the analyst.

A responsible analyst can submit a completion ticket before the due-date:

    {
      analysis: "exome_pipe",
      name: "dt_lung_tumor_pair"
    }

If the loader check passes or a completion ticket is received (or the job is canceled by the user), Polyphemus archives the analysis request.

In the final event, the analyst is responsible for acquiring data and uploading it, presumably directly from/to magma, perhaps to a given data loader. Polyphemus merely catalogues and suggests (and sometimes deletes), he does not insert new data into magma.

None of this precludes loading data into magma 'by hand', but Polyphemus provides a simple dispatch interface to allow a wide variety of users to present and consume work (computational analysis).
