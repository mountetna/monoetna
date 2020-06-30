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

### Managing Biological Data Files

Using Timur, you can view and edit the biological assay files that relate to each Magma record. For example, FastQ files, Flow Cytometry files, biopsy images, etc. These types of attributes will appear with one of three states, depending on the value assigned.

[[/assets/images/timur/sample-view-only.png|Viewing sample attributes on Timur]]

When editing the sample data, users will see three buttons to set the state of the attribute.

[[/assets/images/timur/sample-edit-view.png|Editing sample attributes on Timur reveals three buttons for file and image attributes]]

#### No File

[[/assets/images/timur/no-file-present.png|No File text when file has not been uploaded]]

When no file has been uploaded, the attribute shows a red `No File` value. This means that a file has not been uploaded for the attribute, and that one is expected to be present.

[[/assets/images/timur/delete-file-button-callout.png|Delete button circled in the UI]]

When editing, this action deletes any currently assigned value to the attribute.

#### Blank File

[[/assets/images/timur/closeup-three-file-buttons.png|Blank File text when attribute has been marked as blank]]

When the attribute has explicitly been marked as "Blank", the attribute shows a gray `Blank File` value. This means that the system does not expect a file for the attribute, and it has been deliberately marked as such. Note that this differs from the `No File` value.

[[/assets/images/timur/closeup-three-file-buttons-callout.png|Blank file button circled in the UI]]

When editing, this action marks the file or image attribute as deliberately blank.

#### Link to Metis File

[[/assets/images/timur/link-to-view-file.png|Existing image attribute value links to the file on Metis]]

If a file or image has been assigned to the attribute, a system-generated filename appears as the value and acts as a link. Clicking the link will open the file on Metis.

[[/assets/images/timur/metis-path-button.png|Editing a Metis link opens up a new input text field]]

In edit mode, clicking the Metis link button opens up a new input text field for a Metis path. You can grab a Metis path for an existing file by navigating to Metis, and clicking the right-hand menu for a file. Select the "Copy metis path" option, which will put the path into your clipboard.

[[/assets/images/metis/file-dropdown-copy-path.png|Copy metis path option in dropdown menu]]

Paste this value into the Timur input field, and when you save the record the value will be updated.

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
