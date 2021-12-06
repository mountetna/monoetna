---
layout: default
group: app
---

# Vulcan
{:.no_toc}

## Calculation Service
{:.no_toc}

Vulcan provides tools for doing analysis of data in Mount Etna.

* TOC
{:toc}

### Workflows

A Vulcan workflow is defined in the Common Workflow Language (CWL,
https://www.commonwl.org/user_guide) and runs steps which run Python code or
mediate user interactions.  Users may interact with workflows in the Vulcan UI,
setting input parameters and queueing the workflow to be run by Vulcan.
Subsequently the workflow is executed by Vulcan in a sandboxed Python
environment called Archimedes.

A workflow consists of a set of steps, each with defined inputs and outputs.
The workflow begins with **primary inputs**, and each subsequent step depends
on either primary inputs or the output of a previous step (producing a directed
acyclic graph or DAG).

Here is an example workflow, written in CWL.

     inputs:
       anotherInt:
         type: int
       someInt:
         type: int
     
     steps:
       calcSum:
         run: scripts/add.cwl
         in:
           a: someInt
           b: anotherInt
         out: [sum]
       pickANum:
         run: ui-queries/int.cwl
         in:
           num: calcSum/sum
         out: [num]
       showSum:
         run: ui-outputs/raw.cwl
         in:
           a: pickANum/num
         out: []

The first step in the example makes the two primary inputs `someInt` and
`anotherInt` available to a python script as inputs `a` and `b`. Here is the
corresponding Python script `add.py`:

    from archimedes.functions.dataflow import output_path, input_path

    a = int(open(input_path('a'), 'r').read())
    b = int(open(input_path('b'), 'r').read())

    with open(output_path('sum'), 'w') as output_file:
        output_file.write(str(a + b))

When the step writes output to the appropriate file location, it becomes
available to subsequent (dependent) steps in the workflow. It also becomes part
of Vulcan's cache.

The second step is a UI query step, which uses output from previous steps to
request input from the user using one of Vulcan's query widgets (e.g., a text
box with a label returning an integer), allowing the user to make a response
and write it into the cache for subsequent steps.  There is no associated
Python code; the Vulcan orchestrator will pass the appropriate instruction to
the UI for `ui-queries`.

The last step in the example workflow generates output instructions for the
Vulcan UI, which can be any of a number of output artifacts: an image, a file
download, a table, an interactive plot, etc.

### Vulcan Cache

The cache is an ephemeral data store that holds all of the output produced by
running workflows. Each step in the workflow, when it completes, stores data in
a cache location.  On subsequent runs, if there is already a cached output for
the step, it will be used instead of rerunning the step.  By making use of
cached intermediates, Vulcan workflows can minimize the amount of calculation
required by iteration.

If the inputs or the script for a step change, the location of the step cache
will also change, and Vulcan will look for output in the new location. Note
that the "new" location may also be a previously cached state (e.g., if you
toggle back and forth between two input values). In addition each workflow
session contains a unique key, which if reset will change the cache locations.

### Running a Workflow

Workflows can be run in the Vulcan UI. To run a workflow, navigate to your
project page from the Vulcan root page. Select one of the available workflows
and the Vulcan Workflow UI will appear with a new workflow session.

The interface is divided into three columns. On the left, the user may
configure primary inputs, or, once the workflow is running, UI query inputs
created by the workflow. Outputs, such as plots or download links, appear in
the second column. The final column shows the steps in the workflow in the
order they will run (approximately, since the workflow may branch).

The user may change inputs as they see fit; changes will not be persisted to
the session until the user saves their edits to each required input using the
"Commit" button. Once the session is updated with new inputs, Vulcan will note
that some steps in the workflow can be run to produce new outputs and will
enable the "Run" button.

Upon running, the workflow will run until it completes, encounters an error, or
requires further input from the user to continue.
