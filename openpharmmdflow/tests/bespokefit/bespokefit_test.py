from openff.utilities.testing import skip_if_missing


@skip_if_missing("openff.bespokefit")
def test_bespokefit():
    from openff.bespokefit.executor import BespokeExecutor
    from openff.bespokefit.executor import BespokeWorkerConfig
    from openff.bespokefit.executor.client import BespokeFitClient
    from openff.bespokefit.executor.client import Settings
    from openff.bespokefit.workflows import BespokeWorkflowFactory
    from openff.qcsubmit.common_structures import QCSpec
    from openff.toolkit.topology import Molecule

    factory = BespokeWorkflowFactory(
        # Define the starting force field that will be augmented with bespoke
        # parameters.
        initial_force_field="openff-2.2.1.offxml",
        # Change the level of theory that the reference QC data is generated at
        default_qc_specs=[
            QCSpec(
                method="gfn2xtb",
                basis=None,
                program="xtb",
                spec_name="xtb",
                spec_description="gfn2xtb",
            )
        ],
    )

    input_molecule = Molecule.from_smiles("C(C(=O)O)N")  # Glycine

    workflow_schema = factory.optimization_schema_from_molecule(molecule=input_molecule)

    # create a client to interface with the executor
    settings = Settings()
    client = BespokeFitClient(settings=settings)

    with BespokeExecutor(
        n_fragmenter_workers=1,
        n_optimizer_workers=1,
        n_qc_compute_workers=2,
        qc_compute_worker_config=BespokeWorkerConfig(n_cores=1),
    ) as executor:
        # Submit our workflow to the executor
        task_id = client.submit_optimization(input_schema=workflow_schema)
        # Wait until the executor is done
        output = client.wait_until_complete(task_id)

    if output.status == "success":
        # Save the resulting force field to an OFFXML file
        output.bespoke_force_field.to_file("output-ff.offxml")
    elif output.status == "errored":
        # OR the print the error message if unsuccessful
        print(output.error)

    factory.to_file("workflow-factory.yaml")
    factory = BespokeWorkflowFactory.from_file("workflow-factory.yaml")
