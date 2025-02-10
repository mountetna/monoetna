## Vulcan

Vulcan is the Data Library app powering visualization, processing, and exploration of data.
To achieve this goal, it relies on pre-defined 'workflows' and allows users to create 'workspaces' where a particular version of the workflow can be run.
While users interact with Vulcan from their browser, actual calculation steps of the workflow are run on a remote computation server, and inside of the file directory associated with the given workspace.

## The Vulcan Interface:

The interface is set up in 3 sections:

- User Inputs on the left:
  - The left side is where any interface for setting parameters of the workflow will show up.
  - "Config Params": Most workflows will start with a set, or sets, of Config Params that establish initial settings of the workflow. These might, for example, establish what data to target.
  - Additional parameterization interfaces might show up later on, depending on the design of the workflow.  Generally, when a new interface shows up, you will need to fill it out in order to continue forwards with the workflow.
  - After setting these inputs up as desired, click the `Confirm` button to have your parameterzations sent over to the calculation server  
  - Help Text: Often a circled `?` will appear to the right of a given input. Hovering your cursor over this `?` will display text that the workflow author added in order to describe a particular paramater. The help text may also guide use of the associated user-interface inputs can take many forms and some can be quite complicated.
  - Some example User-Interfaces you migth encounter:
    - simple textbox
    - slider for selecting a number
    - dropdown for selecting one of a set of options
    - multiple dropdowns that power delving into metis folders and selecting a file
    - multiple dropdowns, checkboxes, sliders, and textboxes that together establish the entire parameterization of a visualization
- Outputs in the middle:
  - The center section is where any outputs will appear after successful generation.
  - Generally, you will have to run through multiple steps of a workflow before any outputs will be generated.
  - Outputs can include:
    - visualizations
    - text outputs
    - file download links
- Progress Tracker on the right:
  - The right side shows the set of individual steps of the workflow, each with icons indicating their status. iconsection is where the user can track underlying steps of the pipeline which have been run or are currently running.

## Moving forwards through a workflow

To move forwards through a workflow:

1. Fill out any inputs on the left and use the `Confirm`-button to record them when you are done.
2. Use the `Run` button to trigger any runnable steps of the workflow to be queued, and then run, on the computation server.
3. Wait for completion of steps.  The page should automatically update as steps' change status.  Additional user inputs of outputs may appear as steps complete.
4. Repeat steps 1-3 when new inputs appear, or if you wish to make changes to old parameterizations, until all work is complete.

## Additional Features

To be filled out more... tags, workspace names, revisions, etc.
