# ADR-0004: Gnomon Features

## Status

Proposed

## Context

Assigning and managing names to data entities is a vital function of the
library. To facilitate the process of data entering into the library, it is
necessary to have a well-defined identifier scheme, but also well-established
workflows allowing the generation of new identifiers. To this end we have
created Gnomon, an application which allows us to define identifier grammars,
decompose and compose identifiers, and calculate parent identifiers based on
the identifier decomposition.

However, a number of identifier-related workflows should be developed to allow
Gnomon to be used seamlessly and without significant admin interventions
required to resolve complications:

- Attaching identifier-related metadata to identifiers upon creation, so that this
  data may be easily propagated along with subsequent data into the library
- Deprecating identifiers that are known to be incorrect; data associated with
  these identifiers should be removed from the library.
- Renaming identifiers in all cases where names are incorrect and need to be
  altered; data associated with these identifiers should be put in the correct
  place.
- Alias identifiers generated externally into the data library easily

## Decision

### Attaching metadata

A new section of the Gnomon grammar will define metadata associated with rules.
The rule format will be amended to include an array of required attribute_names
associated with each rule. These attribute names are now required metadata that
must be included when an identifier is created; the endpoint for creating a
name will enforce this requirement. They need not be defined in the magma model
matching the rule at the time of grammar creation, but they must be defined at
the time of identifier creation, or creation will fail. Table or collection
attributes cannot be attached this way.

The Gnomon UI will show, for each identifier, inputs for each required
attribute, with current values if they already exist. The UI will prevent
identifier creation unless each of the required attributes are filled in.

Both data and link attributes may be required. Since link attributes must point
to an existing record, Gnomon will validate that the new identifier exists
before setting them. In some cases, you may create a reciprocal requirement
(that is, you may require the attribute rna_seq.rna_seq_plate and the attribute
rna_seq_plate.rna_seq); the UI will reconcile this reciprocity.

Once Gnomon validates the identifier update, the metadata is deposited into
magma; Gnomon will NOT attach the record if it is not already attached. The
data will thus not become visible until other generated data associated with
the identifier is put into Magma.

After identifier creation, when viewing identifiers as a table, the Gnomon UI
will also allow required metadata to be shown.

### Aliasing identifiers

Aliasing may be considered a special case of attaching metadata, but is
distinct enough that it deserves separate mention. Many projects begin
externally and will come with their own identifiers, either in bulk or as a
stream, or both. Gnomon can support more easily moving data from legacy or
external study identifiers by allowing the alias to be recorded alongside the
identifier creation.

This may be achieved by using the metadata attachment workflow along with an
attribute to hold the external study identifier (or several). If so desired,
these attributes may be restricted to allow them to be recorded without being
visible to researchers. In this case, the user creating identifiers must have
restricted permission to the project.

### Deprecating identifiers

Gnomon will allow identifiers to be deprecated. The identifier record is not
removed, but simply marked as deprecated. For most purposes, Gnomon will
consider this to be a "created" identifier.

When an identifier is marked as deprecated, the record associated with it is
removed from Magma, along with linked files from Metis (data blocks are
destroyed). Subsequent updates to this identifier in Magma will not result in
data being attached; the data will be ignored (non-files) or destroyed from
Metis (files) instead. Therefore deprecation should be understood as a serious
operation that will result in the removal of data.

If it is discovered that an identifier was deprecated in error, or if a new
requirement appears to attach data at that identifier, the deprecation may be
reversed. However, this will not undo the removal of any data that may have
resulted from deprecation.

### Renaming identifiers

Gnomon allows an identifier to be renamed to a new identifier. It is presumed
that there is no data at the new identifier; if Gnomon discovers attached data
at the new identifier, it will not allow renaming. Upon renaming, attached data
and its children will be renamed according to the new identifier and its
descendants.

There are a few scenarios where renaming may occur:
1. Two different samples were given the same name; data was generated for each
   under this name. In this case, the data will have to be manually
   disambiguated and renamed.
2. A sample was given the incorrect name; data was generated under this
   incorrect name. In this case, the (linked) data will be automatically moved
   by Gnomon - unattached files in Metis will not be touched and will remain
   under the old name (a separate workflow may warn of these).
3. One sample was given two different names. Data was generated twice for the
   same sample under each name. One name should be deprecated; data can be
   salvaged manually as a replicate if desired.
4. Two samples were swapped (each given the name of the other) and data
   generated for each under the wrong name. In this case, to avoid further
   confusion in resolution, each swapped name should be given a new identifier.

## Consequences

These changes will allow Gnomon to be more useful in managing the creation and
destruction of data associated with study entities by facilitating identifier
creation.

1. Attaching metadata allows several different kinds of crucial information to
   be associated with samples upon identifier creation. The most prominent
   examples are external study identifiers (aliases) and pool identifiers. Both
   of these are vital to record in association with identifiers; being able to
   attach these at the moment the identifier is defined will greatly increase
   data integrity and ease metadata transfer.
2. Allowing identifier deprecation enables a basic editing function that makes
   current identifier creation cumbersome (it is impossible to undo mistakes).
   In addition, it allows a convenient way to expunge data upon the removal of
   patient consent, a common and vital use case.
3. Renaming identifiers allows us to recover from errors and make sure that
   data attached under incorrect identifiers is correctly propagated to its new
   place.

Overall we hope that these features will make Gnomon an attractive and simple
way to manage the identifier lifecycle and promote the easy attachment of data.
