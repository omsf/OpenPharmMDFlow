"""
Wrapper around bespokefit
"""

from openff.bespokefit.executor import BespokeExecutor
from openff.bespokefit.executor import BespokeWorkerConfig
from openff.bespokefit.executor.client import BespokeFitClient
from openff.bespokefit.executor.client import Settings
from openff.bespokefit.workflows import BespokeWorkflowFactory
from openff.qcsubmit.common_structures import QCSpec


def build_bespoke_workflow_factory(bespoke_workflow_factory_config):
    factory = BespokeWorkflowFactory(
        # Define the starting force field that will be augmented with bespoke
        # parameters.
        initial_force_field=bespoke_workflow_factory_config.initial_force_field,
        # Change the level of theory that the reference QC data is generated at
        default_qc_specs=[
            QCSpec(
                method=bespoke_workflow_factory_config.qc_method,
                basis=bespoke_workflow_factory_config.qc_basis,
                program=bespoke_workflow_factory_config.qc_program,
                spec_name=bespoke_workflow_factory_config.qc_spec_name,
                spec_description=bespoke_workflow_factory_config.qc_spec_description,
            )
        ],
    )
    return factory


def run_bespokefit(bespokefit_config, molecule, bespoke_workflow_factory):
    # reduce nesting
    bespoke_executor_config = bespokefit_config.bespoke_executor_config

    # TODO: this takes in a list of mols, but I am not sure how to correctly
    # pull out the ff
    workflow_schemas = bespoke_workflow_factory.optimization_schemas_from_molecules(
        molecules=molecule, processors=bespoke_executor_config.n_fragment_workers
    )
    # create a client to interface with the executor
    settings = Settings()
    client = BespokeFitClient(settings=settings)

    with BespokeExecutor(
        n_fragmenter_workers=bespoke_executor_config.n_fragmenter_workers,
        n_optimizer_workers=bespoke_executor_config.n_optimizer_workers,
        n_qc_compute_workers=bespoke_executor_config.n_qc_compute_workers,
        qc_compute_worker_config=BespokeWorkerConfig(
            n_cores=bespoke_executor_config.n_bespoke_workers
        ),
    ) as executor:
        # Submit our workflow to the executor
        for workflow_schema in workflow_schemas:
            task_id = client.submit_optimization(input_schema=workflow_schema)
        # Wait until the executor is done
        output = client.wait_until_complete(task_id)

    if output.status == "success":
        # Save the resulting force field to an OFFXML file
        if bespokefit_config.save_bespoke_ff:
            output.bespoke_force_field.to_file("output-ff.offxml")
        return output.bespoke_force_field
    elif output.status == "errored":
        # OR the print the error message if unsuccessful
        print(output.error)
