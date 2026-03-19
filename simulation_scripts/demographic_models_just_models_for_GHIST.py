'''
Add these models to:
stdpopsim/stdpopsim/catalog/DroMel/demographic_models.py
'''

def _afr_growth_testing_10kPop():
    id = "AfricanGrowth_testing_10kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    # mutation_rate = 8.4e-8
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    N_A = 11527
    # # Times are provided in 4N_ref generations, so we convert into generations.
    # # generation_time = 10 / year
    t_1 = 2013.53

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_testing_10kPop())



def _afr_growth_final_10kPop():
    id = "AfricanGrowth_final_10kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    N_A = 9731
    # Times are provided in 4N_ref generations, so we convert into generations.
    # generation_time = 10 / year
    t_1 = 1532.3345 # Generation time

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_final_10kPop())




def _afr_growth_testing_20kPop():
    id = "AfricanGrowth_testing_20kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    # N_ref = 12503
    # t_1_coal = 0.5 / 2
    # t_2_coal = 1 / 2
    # # estimates from the ANN
    # N_R = 15442
    # N_B = 114530
    N_A = 21527
    N_G = 31234.2124
    # # Times are provided in 4N_ref generations, so we convert into generations.
    # # generation_time = 10 / year
    # t_1 = t_1_coal * 4 * N_ref
    # t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    # print(t_1, t_2)

    t_1 = 3013.53
    # t_1_years = t_1 * 0.1
    t_dadi = t_1/(2*N_A)
    nu = N_G/N_A
    tnu = t_dadi/nu

    print(f"dadi nu, one-shot challenge: {nu}")
    print(f"dadi T, one-shot challenge: {t_dadi}")
    print(f"dadi T/nu, one-shot challenge: {tnu}")

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            # # Size change at bottleneck (back in time; BIT)
            # msprime.PopulationParametersChange(
            #     time=t_1, initial_size=N_B, population_id=0
            # ),
            # # Size change at recovery (BIT)
            # msprime.PopulationParametersChange(
            #     time=t_2, initial_size=N_A, population_id=0
            # ),
            # Growth from N_R to N_A
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_testing_20kPop())



def _afr_growth_final_20kPop():
    id = "AfricanGrowth_final_20kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    # # N_ref = 15731
    # t_1_coal = 0.5 / 2
    # t_2_coal = 1 / 2
    # estimates from the ANN
    # N_R = 15442
    # N_B = 114530
    N_A = 19731
    N_G = 34548.86
    # Times are provided in 4N_ref generations, so we convert into generations.
    # generation_time = 10 / year
    # t_1 = t_1_coal * 4 * N_ref
    # t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    # print(t_1, t_2)

    t_1 = 3532.3345 # Generation time
    # t_1_years = t_1 * 0.1
    t_dadi = t_1/(2*N_A)
    nu = N_G/N_A
    tnu = t_dadi/nu

    print(f"dadi nu, one-shot challenge: {nu}")
    print(f"dadi T, one-shot challenge: {t_dadi}")
    print(f"dadi T/nu, one-shot challenge: {tnu}")

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            # # Size change at bottleneck (back in time; BIT)
            # msprime.PopulationParametersChange(
            #     time=t_1, initial_size=N_B, population_id=0
            # ),
            # # Size change at recovery (BIT)
            # msprime.PopulationParametersChange(
            #     time=t_2, initial_size=N_A, population_id=0
            # ),
            # Growth from N_R to N_A
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_final_20kPop())




def _afr_growth_testing_50kPop():
    id = "AfricanGrowth_testing_50kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    # N_ref = 12503
    # t_1_coal = 0.5 / 2
    # t_2_coal = 1 / 2
    # # estimates from the ANN
    # N_R = 15442
    # N_B = 114530
    N_A = 41527
    N_G = 51234.2124
    # # Times are provided in 4N_ref generations, so we convert into generations.
    # # generation_time = 10 / year
    # t_1 = t_1_coal * 4 * N_ref
    # t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    # print(t_1, t_2)

    t_1 = 3013.53
    # t_1_years = t_1 * 0.1
    t_dadi = t_1/(2*N_A)
    nu = N_G/N_A
    tnu = t_dadi/nu

    print(f"dadi nu, one-shot challenge: {nu}")
    print(f"dadi T, one-shot challenge: {t_dadi}")
    print(f"dadi T/nu, one-shot challenge: {tnu}")

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            # # Size change at bottleneck (back in time; BIT)
            # msprime.PopulationParametersChange(
            #     time=t_1, initial_size=N_B, population_id=0
            # ),
            # # Size change at recovery (BIT)
            # msprime.PopulationParametersChange(
            #     time=t_2, initial_size=N_A, population_id=0
            # ),
            # Growth from N_R to N_A
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_testing_50kPop())



def _afr_growth_final_50kPop():
    id = "AfricanGrowth_final_50kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    # # N_ref = 15731
    # t_1_coal = 0.5 / 2
    # t_2_coal = 1 / 2
    # estimates from the ANN
    # N_R = 15442
    # N_B = 114530
    N_A = 49731
    N_G = 54548.86
    # Times are provided in 4N_ref generations, so we convert into generations.
    # generation_time = 10 / year
    # t_1 = t_1_coal * 4 * N_ref
    # t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    # print(t_1, t_2)

    t_1 = 3532.3345 # Generation time
    # t_1_years = t_1 * 0.1
    t_dadi = t_1/(2*N_A)
    nu = N_G/N_A
    tnu = t_dadi/nu

    print(f"dadi nu, one-shot challenge: {nu}")
    print(f"dadi T, one-shot challenge: {t_dadi}")
    print(f"dadi T/nu, one-shot challenge: {tnu}")

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            # # Size change at bottleneck (back in time; BIT)
            # msprime.PopulationParametersChange(
            #     time=t_1, initial_size=N_B, population_id=0
            # ),
            # # Size change at recovery (BIT)
            # msprime.PopulationParametersChange(
            #     time=t_2, initial_size=N_A, population_id=0
            # ),
            # Growth from N_R to N_A
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_final_50kPop())




def _afr_growth_testing_1kPop():
    id = "AfricanGrowth_testing_1kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    # N_ref = 12503
    # t_1_coal = 0.5 / 2
    # t_2_coal = 1 / 2
    # # estimates from the ANN
    # N_R = 15442
    # N_B = 114530
    N_A = 1427
    N_G = 2324
    # # Times are provided in 4N_ref generations, so we convert into generations.
    # # generation_time = 10 / year
    # t_1 = t_1_coal * 4 * N_ref
    # t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    # print(t_1, t_2)

    t_1 = 2013.53
    t_dadi = t_1/(2*N_A)
    nu = N_G/N_A
    tnu = t_dadi/nu

    print(f"dadi nu, one-shot challenge: {nu}")
    print(f"dadi T, one-shot challenge: {t_dadi}")
    print(f"dadi T/nu, one-shot challenge: {tnu}")

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            # # Size change at bottleneck (back in time; BIT)
            # msprime.PopulationParametersChange(
            #     time=t_1, initial_size=N_B, population_id=0
            # ),
            # # Size change at recovery (BIT)
            # msprime.PopulationParametersChange(
            #     time=t_2, initial_size=N_A, population_id=0
            # ),
            # Growth from N_R to N_A
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_testing_1kPop())



def _afr_growth_final_1kPop():
    id = "AfricanGrowth_final_1kPop"
    description = "Customized African3Epoch_1S16 to a simple growth model"
    long_description = """
        A simple growth model for a single African Drosophila Melanogaster population.
    """
    populations = [_afr_population]
    citations = [
        stdpopsim.Citation(
            author="NaN",
            year=2069,
            doi="https://doi.org/10.1371/journal.pcbi.1004845",
            reasons={stdpopsim.CiteReason.DEM_MODEL},
        )
    ]
    generation_time = _species.generation_time
    mutation_rate = 3.36e-8

    # # Parameter values from "Simulating Data" section
    # # these are assumptions, not estimates
    # # N_ref = 15731
    # t_1_coal = 0.5 / 2
    # t_2_coal = 1 / 2
    # estimates from the ANN
    # N_R = 15442
    # N_B = 114530
    N_A = 1571
    N_G = 4748
    # Times are provided in 4N_ref generations, so we convert into generations.
    # generation_time = 10 / year
    # t_1 = t_1_coal * 4 * N_ref
    # t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

    # print(t_1, t_2)

    t_1 = 1532.3345

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_G, metadata=populations[0].asdict()
            )
        ],
        demographic_events=[
            # # Size change at bottleneck (back in time; BIT)
            # msprime.PopulationParametersChange(
            #     time=t_1, initial_size=N_B, population_id=0
            # ),
            # # Size change at recovery (BIT)
            # msprime.PopulationParametersChange(
            #     time=t_2, initial_size=N_A, population_id=0
            # ),
            # Growth from N_R to N_A
            msprime.PopulationParametersChange(
                time=t_1, initial_size=N_A, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_afr_growth_final_1kPop())



# def _afr_growth():
#     id = "AfricanGrowth"
#     description = "Customized African3Epoch_1S16 to a simple growth model"
#     long_description = """
#         A simple growth model for a single African Drosophila Melanogaster population.
#     """
#     populations = [_afr_population]
#     citations = [
#         stdpopsim.Citation(
#             author="NaN",
#             year=2069,
#             doi="https://doi.org/10.1371/journal.pcbi.1004845",
#             reasons={stdpopsim.CiteReason.DEM_MODEL},
#         )
#     ]
#     generation_time = _species.generation_time
#     mutation_rate = 8.4e-9

#     # Parameter values from "Simulating Data" section
#     # these are assumptions, not estimates
#     N_ref = 12503
#     t_1_coal = 0.5 / 2
#     t_2_coal = 1 / 2
#     # estimates from the ANN
#     # N_R = 15442
#     # N_B = 114530
#     N_A = 14527
#     N_G = 18234
#     # Times are provided in 4N_ref generations, so we convert into generations.
#     # generation_time = 10 / year
#     # t_1 = t_1_coal * 4 * N_ref
#     # t_2 = (t_1_coal + t_2_coal) * 4 * N_ref

#     # print(t_1, t_2)

#     t_1 = 12013.5

#     return stdpopsim.DemographicModel(
#         id=id,
#         description=description,
#         long_description=long_description,
#         populations=populations,
#         citations=citations,
#         generation_time=generation_time,
#         mutation_rate=mutation_rate,
#         population_configurations=[
#             msprime.PopulationConfiguration(
#                 initial_size=N_G, metadata=populations[0].asdict()
#             )
#         ],
#         demographic_events=[
#             # # Size change at bottleneck (back in time; BIT)
#             # msprime.PopulationParametersChange(
#             #     time=t_1, initial_size=N_B, population_id=0
#             # ),
#             # # Size change at recovery (BIT)
#             # msprime.PopulationParametersChange(
#             #     time=t_2, initial_size=N_A, population_id=0
#             # ),
#             # Growth from N_R to N_A
#             msprime.PopulationParametersChange(
#                 time=t_1, initial_size=N_A, population_id=0
#             ),
#         ],
#     )
# _species.add_demographic_model(_afr_growth())


def _ooa_2():
    id = "OutOfAfrica_2L06"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.45e-9  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_A0 = 8.603e06
    t_A0 = 600000  # assuming 10 generations / year
    N_A1 = N_A0 / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_E0 = 1.075e06
    N_E1 = 2200
    t_AE = 58000  # generations
    t_E1 = t_AE - 3400

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_A0, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_E0, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            # Size change at Euro bottleneck
            msprime.PopulationParametersChange(
                time=t_E1, initial_size=N_E1, population_id=1
            ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            # African bottleneck
            msprime.PopulationParametersChange(
                time=t_A0, initial_size=N_A1, population_id=0
            ),
        ],
    )
_species.add_demographic_model(_ooa_2())





# def _split_mig_growth_testing_10kPop():
#     id = "SplitMigGrowth_testing_10kPop"
#     description = "Three epoch model for African and European populations"
#     long_description = """
#         The three epoch (modern, bottleneck, ancestral) model estimated for two
#         Drosophila Melanogaster populations: African (ancestral) and European (derived)
#         from Li and Stephan (2006).
#     """
#     populations = [_afr_population, _eur_population]
#     citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
#     generation_time = _species.generation_time
#     mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

#     # African Parameter values from "Demographic History of the African
#     # Population" section
#     N_AF = 12143.89
#     # t_A0 = 123564.54  # assuming 10 generations / year
#     # N_A1 = N_AF / 5.0

#     # European Parameter values from "Demo History of Euro Population"
#     N_EF = 13356.231
#     # N_E1 = 2200
#     t_AE = 4121.876   # generations split
#     # t_E1 = t_AE - 3400 # Split

#     m = 2.458e-5 # Migration rate

#     return stdpopsim.DemographicModel(
#         id=id,
#         description=description,
#         long_description=long_description,
#         populations=populations,
#         citations=citations,
#         generation_time=generation_time,
#         mutation_rate=mutation_rate,
#         population_configurations=[
#             msprime.PopulationConfiguration(
#                 initial_size=N_AF/2, metadata=populations[0].asdict()
#             ),
#             msprime.PopulationConfiguration(
#                 initial_size=N_EF/2, metadata=populations[1].asdict()
#             ),
#         ],
#         demographic_events=[
#             msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
#             msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
#             # # Size change at Euro bottleneck
#             # msprime.PopulationParametersChange(
#             #     time=t_E1, initial_size=N_E1, population_id=1
#             # ),
#             msprime.PopulationParametersChange(
#                 time=0, initial_size=N_AF, population_id=0
#             ),
#             msprime.PopulationParametersChange(
#                 time=0, initial_size=N_EF, population_id=1
#             ),
#             # Split
#             msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
#             msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
#             msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),

#             # # African bottleneck
#             # msprime.PopulationParametersChange(
#             #     time=t_A0, initial_size=N_A1, population_id=0
#             # ),
#         ],
#     )
# _species.add_demographic_model(_split_mig_growth_testing_10kPop())


# def _split_mig_growth_final_10kPop():
#     id = "SplitMigGrowth_final_10kPop"
#     description = "Three epoch model for African and European populations"
#     long_description = """
#         The three epoch (modern, bottleneck, ancestral) model estimated for two
#         Drosophila Melanogaster populations: African (ancestral) and European (derived)
#         from Li and Stephan (2006).
#     """
#     populations = [_afr_population, _eur_population]
#     citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
#     generation_time = _species.generation_time
#     mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

#     # African Parameter values from "Demographic History of the African
#     # Population" section
#     N_AF = 15976.21
#     # t_A0 = 123564.54  # assuming 10 generations / year
#     # N_A1 = N_AF / 5.0

#     # European Parameter values from "Demo History of Euro Population"
#     N_EF = 11372.345
#     # N_E1 = 2200
#     t_AE = 4121.876  # generations split
#     # t_E1 = t_AE - 3400 # Split

#     m = 5.239e-5 # Migration rate

#     return stdpopsim.DemographicModel(
#         id=id,
#         description=description,
#         long_description=long_description,
#         populations=populations,
#         citations=citations,
#         generation_time=generation_time,
#         mutation_rate=mutation_rate,
#         population_configurations=[
#             msprime.PopulationConfiguration(
#                 initial_size=N_AF, metadata=populations[0].asdict()
#             ),
#             msprime.PopulationConfiguration(
#                 initial_size=N_EF, metadata=populations[1].asdict()
#             ),
#         ],
#         demographic_events=[
#             msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
#             msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
#             # # Size change at Euro bottleneck
#             # msprime.PopulationParametersChange(
#             #     time=t_E1, initial_size=N_E1, population_id=1
#             # ),
#             # Split
#             msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
#             msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
#             msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
#             msprime.PopulationParametersChange(
#                 time=t_AE, initial_size=N_AF/2, population_id=0
#             ),
#             msprime.PopulationParametersChange(
#                 time=t_AE, initial_size=N_EF/2, population_id=1
#             ),
#             # # African bottleneck
#             # msprime.PopulationParametersChange(
#             #     time=t_A0, initial_size=N_A1, population_id=0
#             # ),
#         ],
#     )
# _species.add_demographic_model(_split_mig_growth_final_10kPop())





def _split_mig_testing_10kPop():
    id = "SplitMig_testing_10kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 12143.89
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 13356.231
    # N_E1 = 2200
    t_AE = 4221.876   # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 2.458e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_testing_10kPop())


def _split_mig_final_10kPop():
    id = "SplitMig_final_10kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 15976.21
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 11372.345
    # N_E1 = 2200
    t_AE = 4121.876  # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 5.239e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_final_10kPop())






def _split_mig_testing_20kPop():
    id = "SplitMig_testing_20kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 17143.89
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 25356.231
    # N_E1 = 2200
    t_AE = 4221.876   # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 2.458e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_testing_20kPop())


def _split_mig_final_20kPop():
    id = "SplitMig_final_20kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 15976.21
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 22372.345
    # N_E1 = 2200
    t_AE = 4121.876  # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 5.239e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_final_20kPop())






def _split_mig_testing_30kPop():
    id = "SplitMig_testing_30kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 32143.89
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 35356.231
    # N_E1 = 2200
    t_AE = 4221.876   # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 2.458e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_testing_30kPop())


def _split_mig_final_30kPop():
    id = "SplitMig_final_30kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 28726.89
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 34372.345
    # N_E1 = 2200
    t_AE = 4121.876  # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 5.239e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_final_30kPop())




def _split_mig_testing_50kPop():
    id = "SplitMig_testing_50kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 42143.89
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 55356.231
    # N_E1 = 2200
    t_AE = 4121.876   # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 2.458e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_testing_50kPop())


def _split_mig_final_50kPop():
    id = "SplitMig_final_50kPop"
    description = "Three epoch model for African and European populations"
    long_description = """
        The three epoch (modern, bottleneck, ancestral) model estimated for two
        Drosophila Melanogaster populations: African (ancestral) and European (derived)
        from Li and Stephan (2006).
    """
    populations = [_afr_population, _eur_population]
    citations = [_LiAndStephan.because(stdpopsim.CiteReason.DEM_MODEL)]
    generation_time = _species.generation_time
    # mutation_rate = 3.55e-8
    mutation_rate = 1.015e-8  # using the average mutation rate (see citation Methods)

    # African Parameter values from "Demographic History of the African
    # Population" section
    N_AF = 45976.21
    # t_A0 = 123564.54  # assuming 10 generations / year
    # N_A1 = N_AF / 5.0

    # European Parameter values from "Demo History of Euro Population"
    N_EF = 51372.345
    # N_E1 = 2200
    t_AE = 4121.876  # generations split
    # t_E1 = t_AE - 3400 # Split

    m = 5.239e-5 # Migration rate

    return stdpopsim.DemographicModel(
        id=id,
        description=description,
        long_description=long_description,
        populations=populations,
        citations=citations,
        generation_time=generation_time,
        mutation_rate=mutation_rate,
        population_configurations=[
            msprime.PopulationConfiguration(
                initial_size=N_AF, metadata=populations[0].asdict()
            ),
            msprime.PopulationConfiguration(
                initial_size=N_EF, metadata=populations[1].asdict()
            ),
        ],
        demographic_events=[
            msprime.MigrationRateChange(time=0, rate=m, source=1, dest=0),
            msprime.MigrationRateChange(time=0, rate=m, source=0, dest=1),
            # # Size change at Euro bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_E1, initial_size=N_E1, population_id=1
            # ),
            # Split
            msprime.MassMigration(time=t_AE, source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=1, dest=0),
            msprime.MigrationRateChange(time=t_AE, rate=0, source=0, dest=1),
            # # African bottleneck
            # msprime.PopulationParametersChange(
            #     time=t_A0, initial_size=N_A1, population_id=0
            # ),
        ],
    )
_species.add_demographic_model(_split_mig_final_50kPop())
