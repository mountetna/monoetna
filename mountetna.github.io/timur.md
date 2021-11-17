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

Timur's Browse view provides a way to browse individual records by navigating through the data hierarchy.

Using Timur, you can view and edit any attributes on a given Magma record.

In this example, we are viewing the `Nemean Lion` record.
![Viewing sample attributes on Timur](/assets/images/timur/sample-view-only.png)
{:.image}

After clicking the edit button we can edit the attributes on the record.
![Editing sample attributes on Timur reveals three buttons for file and image attributes](/assets/images/timur/sample-edit-view.png)
{:.image}

Edits to any attributes are not saved to Magma until you hit the save button;
you may cancel all edits (including file uploads) before saving.

### Managing Files

File attributes can be set into one of three states: no file, blank, and a reference to a file in Metis.

#### No File

When no file has been uploaded, the attribute shows a red `No file` value. This means that a file has not been uploaded for the attribute, and that one is expected to be present.
![No File text when file has not been uploaded](/assets/images/timur/no-file-present.png)
{:.image}

When editing, this action deletes any currently assigned value to the attribute.
![Delete button circled in the UI](/assets/images/timur/delete-file-button-callout.png)
{:.image}

#### Blank File

When the file has explicitly been marked as "Blank", the attribute shows a gray `Blank file` value. This means, distinct from 'no file', that we do not *expect* a file here, and it has been deliberately marked as blank.
![Blank File text when attribute has been marked as blank](/assets/images/timur/blank-file-attribute.png)
{:.image}

When editing, this action marks the file attribute as blank.
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

## Map

The Map interface shows how your project's data is organized into models. Clicking on a model brings up its list of attributes, which lets you know what data is being stored in Timur.

## Search

The Search interface allows simple filtering of records by their attribute data and bulk download via TSV.
The filter uses the Magma /retrieve api's [filter syntax]({{ site.baseurl }}{% link magma.md %}#retrieve).

## Query

The Query interface allows complex filtering of records across multiple models and bulk download via TSV.
A help guide that provides detailed explanations of the interface, as well as some project-specific examples, is available below.

{% include querying_guide.html %}

## Manifests (advanced users only)

A Manifest allows data to be collected from multiple sources, usually starting
from Magma, and possibly passing through several calculation services (Rtemis
and Pythia). Manifests are written in a scripting language with R-like syntax
(tentatively dubbed TimurLang). Writing a manifest is necessary for extracting
data for use in plots with Timur.

## Plot (advanced users only)

Plots allow data retrieved using a manifest to be graphed in the browser. Timur features a few basic plot types:

1. XY scatter - plot two numerical variables against each other

2. Cluster/heatmap - Plot N numerical variables against each other.
