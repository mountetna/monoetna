import React from 'react';

import markdown from 'etna-js/utils/markdown';

export default function VulcanHelp() {
  const text = "## Vulcan\n\nVulcan is the Data Library app powering visualization, processing, and exploration of data.\nTo achieve this goal, it relies on pre-defined 'workflows' and allows users to create 'workspaces' where a particular version of the workflow can be run.\nWhile users interact with Vulcan from their browser, actual calculation steps of the workflow are run on a remote computation server, and inside of the file directory associated with the given workspace.\n\n## The Vulcan Interface:\n\nThe interface is set up in 3 sections:\n\n- User Inputs on the left:\n  - The left side is where any interface for setting parameters of the workflow will show up.\n  - \"Config Params\": Most workflows will start with a set, or sets, of Config Params that establish initial settings of the workflow. These might, for example, establish what data to target.\n  - Additional parameterization interfaces might show up later on, depending on the design of the workflow.  Generally, when a new interface shows up, you will need to fill it out in order to continue forwards with the workflow.\n  - After setting these inputs up as desired, click the `Confirm` button to have your parameterzations sent over to the calculation server\n  - Help Text: Often a circled `?` will appear to the right of a given input. Hovering your cursor over this `?` will display text that the workflow author added in order to describe a particular paramater. The help text may also guide use of the associated user-interface inputs can take many forms and some can be quite complicated.\n  - Some example User-Interfaces you migth encounter:\n    - simple textbox\n    - slider for selecting a number\n    - dropdown for selecting one of a set of options\n    - multiple dropdowns that power delving into metis folders and selecting a file\n    - multiple dropdowns, checkboxes, sliders, and textboxes that together establish the entire parameterization of a visualization\n- Outputs in the middle:\n  - The center section is where any outputs will appear after successful generation.\n  - Generally, you will have to run through multiple steps of a workflow before any outputs will be generated.\n  - Outputs can include:\n    - visualizations\n    - text outputs\n    - file download links\n- Progress Tracker on the right:\n  - The right side shows the set of individual steps of the workflow, each with icons indicating their status. iconsection is where the user can track underlying steps of the pipeline which have been run or are currently running.\n\n## Moving forwards through a workflow\n\nTo move forwards through a workflow:\n\n1. Fill out any inputs on the left and use the `Confirm`-button to record them when you are done.\n2. Use the `Run` button to trigger any runnable steps of the workflow to be queued, and then run, on the computation server.\n3. Wait for completion of steps.  The page should automatically update as steps' change status.  Additional user inputs of outputs may appear as steps complete.\n4. Repeat steps 1-3 when new inputs appear, or if you wish to make changes to old parameterizations, until all work is complete.\n\n## Additional Features\n\nTo be filled out more... tags, workspace names, revisions, etc."

  return (
    <div
      className='markdown'
      dangerouslySetInnerHTML={{__html: markdown(text)}}
      style={{
        maxHeight: '65vh',
        overflowY: 'auto'
      }}
    />
  );
}
