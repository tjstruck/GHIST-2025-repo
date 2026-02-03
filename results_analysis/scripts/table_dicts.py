
def demo_dict(chal):
    '''
    Returns values of a dictionary that has challenge 
    as a key and the values are the parameters of the
    challenge.
    
    :param chal: Name of challenge
    '''
    demo_dict = {
                    "admixture": 
                        (
                            ["Modern_pop1_proportion_admixture", "Modern_pop1_generations_since_admixture", "Modern_pop2_proportion_admixture", "Modern_pop2_generations_since_admixture"],
                            {
                                "Modern_pop1_proportion_admixture": "$\\%\\text{admix}_1$", 
                                "Modern_pop1_generations_since_admixture": "$T\\text{admix}_1$", 
                                "Modern_pop2_proportion_admixture": "$\\%\\text{admix}_2$", 
                                "Modern_pop2_generations_since_admixture": "$T\\text{admix}_2$"
                            }
                        ),

                    "bottleneck":
                        (    
                            ["generations", "post_decline_fraction"],
                            {
                                "generations": "$T$",
                                "post_decline_fraction": "$\\nu$"
                            }
                        ),

                    "secondary_contact":
                        (
                            ["mainland_population_size", "island_population_size", "generations_since_split", "generations_since_migration", "migration_rate"],
                            {
                                "mainland_population_size": "$N_\\text{main}$",
                                "island_population_size": "$N_\\text{island}$",
                                "generations_since_split": "$T_\\text{split}$",
                                "generations_since_migration": "$T_\\text{mig}$",
                                "migration_rate": "$m$"
                            }
                        ),

                    "growth":
                        (
                            ["generations", "post_growth_scaling"],
                            {
                                "generations": "$T$",
                                "post_growth_scaling": "$\\nu$"
                            }
                        ),

                    "split_migration":
                        (
                            ["AFR_population_size", "AFR_population_size", "generations_since_split", "migration_rate"],
                            {
                                "AFR_population_size": "$N_\\text{AFR}$",
                                "EUR_population_size": "$N_\\text{EUR}$",
                                "generations_since_split": "$T_\\text{split}$",
                                "migration_rate": "$m$"
                            }
                        )

                    }
    return demo_dict[chal]

def participant_dict(uid):
    '''
    Returns the value of a dictionary that takes uid as a 
    key conneting to a formal name of the participant or the group name.
    
    :param uid: Description
    '''
    participant_dict = {
        'peterlaurin': 'Garud Lab Group',
        'alouette': 'Zhang',
        'avaughn': 'Vaughn',
        'tstentella': 'Stentella and Arndt',
        'huangdaxian': 'Xmon Group',
        'Igelkott': 'Mackintosh',
        'srong': 'Rong',
        'austin.t.daigle': 'Daigle',
        'rgollnisch': 'Gollnisch',
        'lzong': 'Zong',
        'wang0207': 'Wang',
        'solomonsloat': 'Sloat',
        'JiatongLiang': 'Liang'
    }

def method_dict(uid):
    '''
    Returns the value of a dictionary that takes uid as a 
    key connecting to the method used.
    
    :param uid: Description
    '''
    method_dict = {
        'peterlaurin': '',
        'alouette': '',
        'avaughn': '',
        'tstentella': '',
        'huangdaxian': '',
        'Igelkott': '',
        'srong': '',
        'austin.t.daigle': '',
        'rgollnisch': '',
        'lzong': '',
        'wang0207': '',
        'solomonsloat': '',
        'JiatongLiang': ''
    }