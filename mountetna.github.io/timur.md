---
layout: default
group: app
---

# Timur

{:.no_toc}

## Data Browser

{:.no_toc}

Timur is a data browser. It is primarily intended to consume data from Magma, a data warehouse.

There are three ways to interact with Timur:

- TOC
  {:toc}

## Browse

The 'browse' view is intended to allow simple record viewing and editing. Magma
publishes a JSON template describing each model and a JSON document describing
each record; Timur merely renders each document using a view built from this
template. This allows us to browse and edit any record in Magma with a single
generic viewer.

Since Magma publishes certain fixed data types, Timur displays each record by
examining the types of the model and rendering it appropriately (for example an
attribute of type 'Date' would be formatted according to local time when
viewing and would show a Date/Time picker when editing).

Occasionally we wish to add custom attributes to the view for a model, or
change the way a certain attribute is displayed (for example, we may always
want to show the Date in Cairo time). In these cases, Timur will patch the
template and record to describe the new attribute and include additional
required data in the record, before passing it on to the client (web browser),
which renders the template as given. This is especially useful for adding custom
plots to a view; we might add a new attribute of class 'BarPlot' to our model,
and include the appropriate bar plot data in the record, which may be rendered
by a component using d3.js.

### Managing Data Files

<div class="row">
  <div class="column" markdown="1">
  <p>Using Timur, you can view and edit any files that relate to a given Magma record. For example, data files, images, PDFs, etc. These types of attributes will appear with one of three states, depending on the value assigned.</p>

  <p markdown="1">
  In this example, `Stats` and `Avatar` are two file attributes for the `Nemean Lion` record.</p>

  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Viewing sample attributes on Timur](/assets/images/timur/sample-view-only.png)
  {: refdef}
  </div>
</div>

<div class="row">
  <div class="column" markdown="1">
  <p>When editing the sample data, users will see three buttons to set the state of the attribute.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Editing sample attributes on Timur reveals three buttons for file and image attributes](/assets/images/timur/sample-edit-view.png)
  {: refdef}
  </div>
</div>

#### No File

<div class="row">
  <div class="column" markdown="1">
  <p>When no file has been uploaded, the attribute shows a red `No file` value. This means that a file has not been uploaded for the attribute, and that one is expected to be present.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![No File text when file has not been uploaded](/assets/images/timur/no-file-present.png)
  {: refdef}
  </div>
</div>

<div class="row">
  <div class="column" markdown="1">
  <p>When editing, this action deletes any currently assigned value to the attribute.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Delete button circled in the UI](/assets/images/timur/delete-file-button-callout.png)
  {: refdef}
  </div>
</div>

#### Blank File

<div class="row">
  <div class="column" markdown="1">
  <p>When the attribute has explicitly been marked as "Blank", the attribute shows a gray `Blank file` value. This means that the system does not expect a file for the attribute, and it has been deliberately marked as such. Note that this differs from the `No file` value.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Blank File text when attribute has been marked as blank](/assets/images/timur/blank-file-attribute.png)
  {: refdef}
  </div>
</div>

<div class="row">
  <div class="column" markdown="1">
  <p>When editing, this action marks the file attribute as deliberately blank.</p>
  </div>
  <div class="column" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Blank file button circled in the UI](/assets/images/timur/blank-file-attribute-callout.png)
  {: refdef}
  </div>
</div>

#### Link to Metis File

<div class="row">
  <div class="column" markdown="1">
  <p>If a file has been assigned to the attribute, a system-generated filename appears as the value and acts as a link. Clicking the link will open the file on Metis.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Existing image attribute value links to the file on Metis](/assets/images/timur/link-to-view-file.png)
  {: refdef}
  </div>
</div>

<div class="row">
  <div class="column" markdown="1">
  <p>In edit mode, clicking the Metis link button lets you input a path to a file that has already been uploaded to Metis.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Editing a Metis link opens up a modal with a new input text field](/assets/images/timur/metis-path-button.png)
  {: refdef}
  </div>
</div>

<div class="row">
  <div class="column" markdown="1">
  <p>Clicking the Metis link button opens up a modal that includes an input text field for a Metis path.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Modal window with input field for metis path](/assets/images/timur/metis-path-modal.png)
  {: refdef}
  </div>
</div>

<div class="row">
  <div class="column" markdown="1">
  <p>You can grab a Metis path for an existing file by navigating to Metis, and clicking the right-hand menu for a file. Select the "Copy metis path" option, which will put the path into your clipboard.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Copy metis path option in dropdown menu](/assets/images/metis/file-dropdown-copy-path.png)
  {: refdef}
  </div>
</div>

<div class="row">
  <div class="column" markdown="1">
  <p>Paste this value into the Timur input field, and when you save the record the value will be updated.</p>
  </div>
  <div class="column image" markdown="1">
  {:refdef: style="text-align: center;margin: 1rem"}
  ![Modal with metis path pasted into the input text field](/assets/images/timur/metis-path-modal-with-data.png)
  {: refdef}
  </div>
</div>

## Search

The Search interface allows simple filtering of records by their attribute data and bulk download via TSV.
The filter uses the Magma /retrieve api's filter syntax, described in detail here: https://github.com/mountetna/magma/wiki/Retrieve

## Manifests

A Manifest allows data to be collected from multiple sources, usually starting
from Magma, and possibly passing through several calculation services (Rtemis
and Pythia). Manifests are written in a scripting language with R-like syntax
(tentatively dubbed TimurLang). Writing a manifest is necessary for extracting
data for use in plots with Timur.

## Plot

Plots allow data retrieved using a manifest to be graphed in the browser. Timur features a few basic plot types:

1. XY scatter - plot two numerical variables against each other

2. Cluster/heatmap - Plot N numerical variables against each other.
