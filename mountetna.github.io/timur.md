---
layout: default
group: app
---

# Timur
{:.no_toc}

## Data Browser
{:.no_toc}

Timur is a data browser. It is primarily intended to view data in Magma (the data warehouse).

* TOC
{:toc}

## Browse

Timur's Browse view provides a way to browse individual records by navigating up and down the data hierarchy.

### Managing Data Files

Using Timur, you can view and edit any files that relate to a given Magma record. For example, data files, images, PDFs, etc. These types of attributes will appear with one of three states, depending on the value assigned.

In this example, `Stats` and `Avatar` are two file attributes for the `Nemean Lion` record.
![Viewing sample attributes on Timur](/assets/images/timur/sample-view-only.png)
{:.image}


When editing the sample data, users will see three buttons to set the state of the attribute.
![Editing sample attributes on Timur reveals three buttons for file and image attributes](/assets/images/timur/sample-edit-view.png)
{:.image}

#### No File

When no file has been uploaded, the attribute shows a red `No file` value. This means that a file has not been uploaded for the attribute, and that one is expected to be present.
![No File text when file has not been uploaded](/assets/images/timur/no-file-present.png)
{:.image}

When editing, this action deletes any currently assigned value to the attribute.
![Delete button circled in the UI](/assets/images/timur/delete-file-button-callout.png)
{:.image}

#### Blank File

When the attribute has explicitly been marked as "Blank", the attribute shows a gray `Blank file` value. This means that the system does not expect a file for the attribute, and it has been deliberately marked as such. Note that this differs from the `No file` value.
![Blank File text when attribute has been marked as blank](/assets/images/timur/blank-file-attribute.png)
{:.image}

When editing, this action marks the file attribute as deliberately blank.
![Blank file button circled in the UI](/assets/images/timur/blank-file-attribute-callout.png)
{:.image}

#### Link to Metis File

If a file has been assigned to the attribute, a system-generated filename appears as the value and acts as a link. Clicking the link will open the file on Metis.
![Existing image attribute value links to the file on Metis](/assets/images/timur/link-to-view-file.png)
{:.image}

In edit mode, clicking the Metis link button lets you input a path to a file that has already been uploaded to Metis.
![Editing a Metis link opens up a modal with a new input text field](/assets/images/timur/metis-path-button.png)
{:.image}

Clicking the Metis link button opens up a modal that includes an input text field for a Metis path.
![Modal window with input field for metis path](/assets/images/timur/metis-path-modal.png)
{:.image}

You can grab a Metis path for an existing file by navigating to Metis, and clicking the right-hand menu for a file. Select the "Copy metis path" option, which will put the path into your clipboard.
![Copy metis path option in dropdown menu](/assets/images/metis/file-dropdown-copy-path.png)
{:.image}

Paste this value into the Timur input field, and when you save the record the value will be updated.
![Modal with metis path pasted into the input text field](/assets/images/timur/metis-path-modal-with-data.png)
{:.image}

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
