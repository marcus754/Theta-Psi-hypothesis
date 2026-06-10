"""Provisional inference layer for Θ–Ψ.

This package is not the physical kernel. Its role is narrower and provisional:
- load observational inputs;
- request predictions from the kernel through stable interfaces;
- assemble likelihoods, bridge relations and sampling utilities.

Key modules:
- `core_api`: single access layer to the physical kernel;
- `model`: background observables from the FRW module;
- `likelihoods*`: likelihood assembly;
- `run_mcmc`: parameter inference and joint-fit orchestration.
"""
