# Project configuration

`bundles.yml` is the central registry for durable bundle identities, roles,
locations, and parent relationships. Each bundle also has one project-specific
configuration file that records its source boundary and workflow state.

The configuration files describe this repository's inputs and decisions. They
do not replace the scWorkbench operation documentation in the supplied image.
Update a planned bundle only after source inspection establishes the actual
inputs, metadata mappings, feature universes, and supported operations.

Bundle status values used here are:

- `planned`: source selection or bundle creation is still pending;
- `working`: a bundle exists but has not completed recipient validation;
- `accepted`: the recorded bundle state passed the project validation gate.

Large inputs, bundles, reports, and run workspaces belong under the mirrored
storage root recorded in `bundles.yml`, outside Git.

A bundle can declare `runtime.r_library` and `runtime.scworkbench_version` when
it needs a separately built scWorkbench runtime in addition to the supplied
image. The runtime's build evidence must be recorded with the bundle. Parse
processing and figure runners resolve both the active bundle path and runtime
from the registry; their first bundle argument still permits an explicit path.
`validate_registered_bundle()` loads the declared library before opening the
bundle, and refuses an incompatible scWorkbench namespace already loaded in R.
