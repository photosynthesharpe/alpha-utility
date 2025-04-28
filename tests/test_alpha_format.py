"""
Spot checks for alpha_format.py.

Author: Serena G. Lotreck
"""
import pytest
from Bio import SeqIO
import sys

sys.path.append('../formatting')
import alpha_format as af

############################ format_protein_seqs ##############################


@pytest.fixture
def fasta_sequences():
    return SeqIO.index('proteins.fasta', 'fasta')


@pytest.fixture
def anchor_sequences():
    return SeqIO.index('anchor_proteins.fasta', 'fasta')


@pytest.fixture
def multimers():
    return {
        'fake_multimer_base_fasta': {
            'CnPGP': 4,
            'TaPGP': 6
        },
        'fake_multimer_anchor_fasta': {
            'SePGP': 8
        }
    }


@pytest.fixture
def in_complex():
    return {
        'CnPGP': 'fake_multimer_base_fasta',
        'TaPGP': 'fake_multimer_base_fasta',
        'SePGP': 'fake_multimer_anchor_fasta'
    }


@pytest.fixture
def base_name():
    return 'CnPGP'


@pytest.fixture
def anchor_name():
    return 'SePGP'


@pytest.fixture
def without_multimers_result():
    return ([{
        'protein': {
            'id': ['A'],
            'sequence':
            'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
        }
    }], 1)


@pytest.fixture
def with_multimers_base_result():
    return ([{
        'protein': {
            'id': ['A', 'B', 'C', 'D'],
            'sequence':
            'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
        }
    }, {
        'protein': {
            'id': ['E', 'F', 'G', 'H', 'I', 'J'],
            'sequence':
            'MSAAQPLTDADHLIDSVETFIFDCDGVIWKGDKLIDGVPQTLDMLRSKGKRLVFVTNNSTKSRKQYGKKFETLGLNVNEEEIFASSFAAAAYLQSIVFPKDKKVYVIGEDGILKELELAGFQYLGGPEDGDKKIELKPGFLMEHDKDVGAVVVGFDRNFNYYKIQYGTLCIRENPGCLFIATNRDAVTHLTDAQEWAGGGSMVGALCGSTQRQPLVVGKPSTFMMDYLANKFGILKSQICMVGDRLDTDILFGQNGGCKTLLVLSGVTTLSMLQSPDNSIQPDFYTNKISDFLSLKATAL'
        }
    }], 10)


@pytest.fixture
def with_multimers_anchor_result():
    return ([{
        'protein': {
            'id': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
            'sequence':
            'MQAIIFDFDGTLVDSLPTVVAIANAHAPDFGYDPIDERDYAQLRQWSSRTIVRRAGLSPWQQARLLQRVQRQLGDCLPALQLFPGVADLLAQLRSRSLCLGILSSNSRQNIEAFLQRQGLRSLFSVVQAGTPILSKRRALSQLVAREGWQPAAVMYVGDETRDVEAARQVGLIAVAVTWGFNDRQSLVAACPDWLLETPSDLLQAVTQLMRQ'
        }
    }], 8)


def test_format_protein_seqs_without_multimers(base_name, fasta_sequences,
                                               without_multimers_result):

    result = af.format_protein_seqs(base_name, fasta_sequences)

    assert result == without_multimers_result


def test_format_protein_seqs_with_multimers_base(base_name, fasta_sequences,
                                                 multimers, in_complex,
                                                 with_multimers_base_result):

    result = af.format_protein_seqs(base_name,
                                    fasta_sequences,
                                    multi_subunit_proteins=multimers,
                                    in_complex=in_complex)

    assert result == with_multimers_base_result


def test_format_protein_seqs_with_multimers_anchor(
        anchor_name, fasta_sequences, anchor_sequences, multimers, in_complex,
        with_multimers_anchor_result):

    result = af.format_protein_seqs(anchor_name,
                                    fasta_sequences,
                                    anchor_sequences,
                                    multi_subunit_proteins=multimers,
                                    in_complex=in_complex)

    assert result == with_multimers_anchor_result


############################ generate_seq_dicts ###############################


@pytest.fixture
def combo_folding_only():
    return ('P_CnPGP', )


@pytest.fixture
def combo_one_v_many():
    return ('P_CnPGP', 'P_SePGP', 'C_Mg')


@pytest.fixture
def combo_many_v_many():
    return ('P_TaPGP', 'P_CnPGP', 'L_PGA')


@pytest.fixture
def folding_only_result():
    return {
        'name':
        'CnPGP',
        'sequences': [{
            'protein': {
                'id': ['A'],
                'sequence':
                'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
            }
        }],
        'dialect':
        'alphafold3',
        'version':
        2,
        'modelSeeds': [1855]
    }


@pytest.fixture
def one_v_many_result():
    return {
        'name':
        'CnPGP_SePGP_C_Mg',
        'sequences': [{
            'protein': {
                'id': ['A'],
                'sequence':
                'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
            }
        }, {
            'protein': {
                'id': ['B'],
                'sequence':
                'MQAIIFDFDGTLVDSLPTVVAIANAHAPDFGYDPIDERDYAQLRQWSSRTIVRRAGLSPWQQARLLQRVQRQLGDCLPALQLFPGVADLLAQLRSRSLCLGILSSNSRQNIEAFLQRQGLRSLFSVVQAGTPILSKRRALSQLVAREGWQPAAVMYVGDETRDVEAARQVGLIAVAVTWGFNDRQSLVAACPDWLLETPSDLLQAVTQLMRQ'
            }
        }, {
            'ligand': {
                'id': ['C'],
                'ccdCodes': ['Mg']
            }
        }],
        'dialect':
        'alphafold3',
        'version':
        2,
        'modelSeeds': [1855]
    }


@pytest.fixture
def many_v_many_result():
    return {
        'name':
        'TaPGP_CnPGP_PGA',
        'sequences': [{
            'protein': {
                'id': ['A'],
                'sequence':
                'MSAAQPLTDADHLIDSVETFIFDCDGVIWKGDKLIDGVPQTLDMLRSKGKRLVFVTNNSTKSRKQYGKKFETLGLNVNEEEIFASSFAAAAYLQSIVFPKDKKVYVIGEDGILKELELAGFQYLGGPEDGDKKIELKPGFLMEHDKDVGAVVVGFDRNFNYYKIQYGTLCIRENPGCLFIATNRDAVTHLTDAQEWAGGGSMVGALCGSTQRQPLVVGKPSTFMMDYLANKFGILKSQICMVGDRLDTDILFGQNGGCKTLLVLSGVTTLSMLQSPDNSIQPDFYTNKISDFLSLKATAL'
            }
        }, {
            'protein': {
                'id': ['B'],
                'sequence':
                'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
            }
        }, {
            'ligand': {
                'id': ['C'],
                'ccdCodes': ['PGA']
            }
        }],
        'dialect':
        'alphafold3',
        'version':
        2,
        'modelSeeds': [1855]
    }


def test_generate_seq_dicts_folding_only(fasta_sequences, combo_folding_only,
                                         folding_only_result):

    result = af.generate_seq_dicts(combo_folding_only, fasta_sequences)

    assert result == folding_only_result


def test_generate_seq_dicts_one_v_many(fasta_sequences, anchor_sequences,
                                       combo_one_v_many, one_v_many_result):

    result = af.generate_seq_dicts(combo_one_v_many, fasta_sequences,
                                   anchor_sequences)

    assert result == one_v_many_result


def test_generate_seq_dicts_many_v_many(fasta_sequences, combo_many_v_many,
                                        many_v_many_result):

    result = af.generate_seq_dicts(combo_many_v_many, fasta_sequences)

    assert result == many_v_many_result


############################## generateJSONs ##################################
# Note that I'm not testing all the possible input formations here, mainly
# due to the manual labor of formulating the correct answers. Instead, I'm
# focusing on the cases that have some unique feature that I want to test, like
# having no cofactors or the many_v_many scenario that involves needing to use
# the same protein list twice
## TODO probably want to add a test for multiple cofactors if that's a use case
## we expect to see


@pytest.fixture
def ligand_list_multiple():
    return ['PGA', 'GOL']


@pytest.fixture
def cofactor_list_one():
    return ['MG']


@pytest.fixture
def cofactor_list_empty():
    return []


@pytest.fixture
def CnPGP_seq():
    return (
        'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
    )


@pytest.fixture
def TaPGP_seq():
    return (
        'MSAAQPLTDADHLIDSVETFIFDCDGVIWKGDKLIDGVPQTLDMLRSKGKRLVFVTNNSTKSRKQYGKKFETLGLNVNEEEIFASSFAAAAYLQSIVFPKDKKVYVIGEDGILKELELAGFQYLGGPEDGDKKIELKPGFLMEHDKDVGAVVVGFDRNFNYYKIQYGTLCIRENPGCLFIATNRDAVTHLTDAQEWAGGGSMVGALCGSTQRQPLVVGKPSTFMMDYLANKFGILKSQICMVGDRLDTDILFGQNGGCKTLLVLSGVTTLSMLQSPDNSIQPDFYTNKISDFLSLKATAL'
    )


@pytest.fixture
def TfPGP_seq():
    return (
        'MSTAQPLTDADHLIDSVQTFIFDCDGVIWKGDKLIDGVPQTLDMLRSKGKRLVFVTNNSTKSRKQYGKKFETLGLNVNEEEIFASSFAAAAYLQSIDFPKDKKVYVIGEDGILKELELAGFQYLGGPEDGDKKIELKPGFLMEHDKDVGAVVVGFDRYFNYYKIQYGTLCIRENSGCLFIATNRDAVTHLTDAQEWAGGGSMVGALCGSTQRQPLVVGKPSTFMMDYLANKFGILKSQICMVGDRLDTDILFGQNGGCKTLLVLSGVTTLSMLQSPDNSIQPDFYTNKISDFLSLKATAR'
    )


@pytest.fixture
def SePGP_seq():
    return (
        'MQAIIFDFDGTLVDSLPTVVAIANAHAPDFGYDPIDERDYAQLRQWSSRTIVRRAGLSPWQQARLLQRVQRQLGDCLPALQLFPGVADLLAQLRSRSLCLGILSSNSRQNIEAFLQRQGLRSLFSVVQAGTPILSKRRALSQLVAREGWQPAAVMYVGDETRDVEAARQVGLIAVAVTWGFNDRQSLVAACPDWLLETPSDLLQAVTQLMRQ'
    )


@pytest.fixture
def folding_only_ligands_and_cofactor_results(CnPGP_seq, TfPGP_seq, TaPGP_seq):
    return {
        'CnPGP_PGA_MG': {
            'name':
            'CnPGP_PGA_MG',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': CnPGP_seq
                }
            }, {
                'ligand': {
                    'id': ['B'],
                    'ccdCodes': ['PGA']
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['MG']
                }
            }]
        },
        'TaPGP_PGA_MG': {
            'name':
            'TaPGP_PGA_MG',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TaPGP_seq
                }
            }, {
                'ligand': {
                    'id': ['B'],
                    'ccdCodes': ['PGA']
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['MG']
                }
            }]
        },
        'TfPGP_PGA_MG': {
            'name':
            'TfPGP_PGA_MG',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TfPGP_seq
                }
            }, {
                'ligand': {
                    'id': ['B'],
                    'ccdCodes': ['PGA']
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['MG']
                }
            }]
        },
        'CnPGP_GOL_MG': {
            'name':
            'CnPGP_GOL_MG',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': CnPGP_seq
                }
            }, {
                'ligand': {
                    'id': ['B'],
                    'ccdCodes': ['GOL']
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['MG']
                }
            }]
        },
        'TaPGP_GOL_MG': {
            'name':
            'TaPGP_GOL_MG',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TaPGP_seq
                }
            }, {
                'ligand': {
                    'id': ['B'],
                    'ccdCodes': ['GOL']
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['MG']
                }
            }]
        },
        'TfPGP_GOL_MG': {
            'name':
            'TfPGP_PGA_MG',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TfPGP_seq
                }
            }, {
                'ligand': {
                    'id': ['B'],
                    'ccdCodes': ['GOL']
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['MG']
                }
            }]
        }
    }


@pytest.fixture
def one_v_many_ligands_only_results(CnPGP, TaPGP, TfPGP, SePGP):
    return {
        'CnPGP_SePGP_PGA': {
            'name':
            'CnPGP_SePGP_PGA',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': CnPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': SePGP_seq
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['PGA']
                }
            }]
        },
        'TaPGP_SePGP_PGA': {
            'name':
            'TaPGP_SePGP_PGA',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TaPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': SePGP_seq
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['PGA']
                }
            }]
        },
        'TfPGP_SePGP_PGA': {
            'name':
            'TfPGP_SePGP_PGA',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TfPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': SePGP_seq
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['PGA']
                }
            }]
        },
        'CnPGP_SePGP_GOL': {
            'name':
            'CnPGP_SePGP_GOL',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': CnPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': SePGP_seq
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['GOL']
                }
            }]
        },
        'TaPGP_SePGP_GOL': {
            'name':
            'TaPGP_SePGP_GOL',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TaPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': SePGP_seq
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['GOL']
                }
            }]
        },
        'TfPGP_SePGP_GOL': {
            'name':
            'TfPGP_SePGP_GOL',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TfPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': SePGP_seq
                }
            }, {
                'ligand': {
                    'id': ['C'],
                    'ccdCodes': ['GOL']
                }
            }]
        }
    }


@pytest.fixture
def many_v_many_no_ligand_no_cofactor_results(CnPGP_seq, TaPGP_seq, TfPGP_seq):
    return {
        'CnPGP_TaPGP': {
            'name':
            'CnPGP_TaPGP',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': CnPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': TaPGP_seq
                }
            }]
        },
        'CnPGP_TfPGP': {
            'name':
            'CnPGP_TfPGP',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': CnPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': TfPGP_seq
                }
            }]
        },
        'TaPGP_TfPGP': {
            'name':
            'TaPGP_TfPGP',
            'sequences': [{
                'protein': {
                    'id': ['A'],
                    'sequence': TaPGP_seq
                }
            }, {
                'protein': {
                    'id': ['B'],
                    'sequence': TfPGP_seq
                }
            }]
        }
    }


def test_generateJSONs_folding_only_ligand_and_cofactor(
        fasta_sequences, ligand_list_multiple, cofactor_list_one,
        folding_only_ligands_and_cofactor_results):

    result = af.generateJSONs(fasta_sequences,
                              ligand_list=ligand_list_multiple,
                              cofactor_list=cofactor_list_one)

    assert result == folding_only_ligands_and_cofactor_results


def test_generateJSONs_one_v_many_ligands_only(
        fasta_sequences, ligand_list_multiple, anchor_sequences,
        one_v_many_ligands_only_results):

    result = af.generateJSONs(fasta_sequences,
                              ligand_list=ligand_list_multiple,
                              anchor_sequences=anchor_sequences,
                              protein_comparison_type='one_v_many')

    assert result == one_v_many_ligands_only_results


def test_generateJSONs_one_v_many_ligands_only(
        fasta_sequences, cofactor_list_empty,
        many_v_many_no_ligand_no_cofactor_results):

    result = af.generateJSONs(fasta_sequences,
                              cofactor_list=cofactor_list_empty,
                              protein_comparison_type='many_v_many')

    assert result == many_v_many_no_ligand_no_cofactor_results
