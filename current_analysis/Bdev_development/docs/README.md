# Development documentation

Use this directory for the scWorkbench capability-development record. For each
Bdev requirement, document:

1. the scientific contract and why an existing operation is insufficient;
2. the owning component: core, an existing official pack, a new official pack,
   or approved outside preprocessing;
3. typed inputs, parameters, outputs, missing-value behavior, and provenance;
4. alternative methods as separate qualified operations when their scientific
   interpretations differ;
5. synthetic contract tests and Bdev acceptance tests;
6. the source branch, reviewed commit, package version, and manifest identity;
7. the SIF build in which the operation is qualified.

Core changes must remain generic. Scientific method implementations belong in
operation packs. Any proposed scWorkbench or PyBench change requires a written
motivation showing the missing core capability, why a pack cannot solve it, the
minimal change, and its qualification impact. Stop and obtain explicit user
approval for that specific proposal before creating a branch or changing core
files. The implementation-stage analysis plan is the requirements source for
this work.
