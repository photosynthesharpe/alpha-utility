"""
Spot checks for alpha_format.py.

Author: Serena G. Lotreck
"""
import pytest
from Bio import SeqIO
import sys

sys.path.append('../formatting')
import alpha_format as af

############################ generate_seq_dicts ###############################


@pytest.fixture
def fasta_sequences():
    return SeqIO.index('proteins.fasta', 'fasta')


@pytest.fixture
def anchor_sequences():
    return SeqIO.index('anchor_proteins.fasta', 'fasta')


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
        'P_CnPGP',
        'sequences': [{
            'protein': {
                'id':
                'P',
                'sequence':
                'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
            }
        }]
    }


@pytest.fixture
def one_v_many_result():
    return {
        'name':
        'P_CnPGP_P_SePGP_C_Mg',
        'sequences': [{
            'protein': {
                'id':
                'P',
                'sequence':
                'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
            }
        }, {
            'protein': {
                'id':
                'P',
                'sequence':
                'MQAIIFDFDGTLVDSLPTVVAIANAHAPDFGYDPIDERDYAQLRQWSSRTIVRRAGLSPWQQARLLQRVQRQLGDCLPALQLFPGVADLLAQLRSRSLCLGILSSNSRQNIEAFLQRQGLRSLFSVVQAGTPILSKRRALSQLVAREGWQPAAVMYVGDETRDVEAARQVGLIAVAVTWGFNDRQSLVAACPDWLLETPSDLLQAVTQLMRQ'
            }
        }, {
            'ligand': {
                'id': 'C',
                'ccdCodes': ['Mg']
            }
        }]
    }


@pytest.fixture
def many_v_many_result():
    return {
        'name':
        'P_TaPGP_P_CnPGP_L_PGA',
        'sequences': [{
            'protein': {
                'id':
                'P',
                'sequence':
                'MSAAQPLTDADHLIDSVETFIFDCDGVIWKGDKLIDGVPQTLDMLRSKGKRLVFVTNNSTKSRKQYGKKFETLGLNVNEEEIFASSFAAAAYLQSIVFPKDKKVYVIGEDGILKELELAGFQYLGGPEDGDKKIELKPGFLMEHDKDVGAVVVGFDRNFNYYKIQYGTLCIRENPGCLFIATNRDAVTHLTDAQEWAGGGSMVGALCGSTQRQPLVVGKPSTFMMDYLANKFGILKSQICMVGDRLDTDILFGQNGGCKTLLVLSGVTTLSMLQSPDNSIQPDFYTNKISDFLSLKATAL'
            }
        }, {
            'protein': {
                'id':
                'P',
                'sequence':
                'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
            }
        }, {
            'ligand': {
                'id': 'L',
                'ccdCodes': ['PGA']
            }
        }]
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
def CnPGP():
    return {
        'protein': {
            'id':
            'P',
            'sequence':
            'MGRATVDDKKALVDKVDCFIFDCDGVIWKGDSVIDGVPETLDMLRKLNKRLVFVTNNSTKSRAGYLGKFTSLGLKVKAEEIYSSSYAAAAYLESINFKKKVYVVGEVGIQEELDLKGISHLGGPADADKKVTLKEGVFFGHDHEVGAVVVGFDRNINYHKIQYATLCIRENPGCLFIATNRDAVTHLTEAQEWAGNGSMVGAIIGSTKREPITVGKPNGFMLENIAKSYGLKPEQICMVGDRLDTDIMFGKNGGLTTCLVLSGVTTEEELLSPKNTIAPDFYMNQLSDMLAIQNSVGSYVEA'
        }
    }


@pytest.fixture
def TaPGP():
    return {
        'protein': {
            'id':
            'P',
            'sequence':
            'MSAAQPLTDADHLIDSVETFIFDCDGVIWKGDKLIDGVPQTLDMLRSKGKRLVFVTNNSTKSRKQYGKKFETLGLNVNEEEIFASSFAAAAYLQSIVFPKDKKVYVIGEDGILKELELAGFQYLGGPEDGDKKIELKPGFLMEHDKDVGAVVVGFDRNFNYYKIQYGTLCIRENPGCLFIATNRDAVTHLTDAQEWAGGGSMVGALCGSTQRQPLVVGKPSTFMMDYLANKFGILKSQICMVGDRLDTDILFGQNGGCKTLLVLSGVTTLSMLQSPDNSIQPDFYTNKISDFLSLKATAL'
        }
    }


@pytest.fixture
def TfPGP():
    return {
        'protein': {
            'id':
            'P',
            'sequence':
            'MSTAQPLTDADHLIDSVQTFIFDCDGVIWKGDKLIDGVPQTLDMLRSKGKRLVFVTNNSTKSRKQYGKKFETLGLNVNEEEIFASSFAAAAYLQSIDFPKDKKVYVIGEDGILKELELAGFQYLGGPEDGDKKIELKPGFLMEHDKDVGAVVVGFDRYFNYYKIQYGTLCIRENSGCLFIATNRDAVTHLTDAQEWAGGGSMVGALCGSTQRQPLVVGKPSTFMMDYLANKFGILKSQICMVGDRLDTDILFGQNGGCKTLLVLSGVTTLSMLQSPDNSIQPDFYTNKISDFLSLKATAR'
        }
    }


@pytest.fixture
def SePGP():
    return {
        'protein': {
            'id':
            'P',
            'sequence':
            'MQAIIFDFDGTLVDSLPTVVAIANAHAPDFGYDPIDERDYAQLRQWSSRTIVRRAGLSPWQQARLLQRVQRQLGDCLPALQLFPGVADLLAQLRSRSLCLGILSSNSRQNIEAFLQRQGLRSLFSVVQAGTPILSKRRALSQLVAREGWQPAAVMYVGDETRDVEAARQVGLIAVAVTWGFNDRQSLVAACPDWLLETPSDLLQAVTQLMRQ'
        }
    }


@pytest.fixture
def pga():
    return {'ligand': {'id': 'L', 'ccdCodes': ['PGA']}}


@pytest.fixture
def gol():
    return {'ligand': {'id': 'L', 'ccdCodes': ['GOL']}}


@pytest.fixture
def mg():
    return {'ligand': {'id': 'C', 'ccdCodes': ['MG']}}


@pytest.fixture
def folding_only_ligands_and_cofactor_results(CnPGP, TfPGP, TaPGP, pga, gol,
                                              mg):
    return {
        'P_CnPGP_L_PGA_C_MG': {
            'name': 'P_CnPGP_L_PGA_C_MG',
            'sequences': [CnPGP, pga, mg]
        },
        'P_TaPGP_L_PGA_C_MG': {
            'name': 'P_TaPGP_L_PGA_C_MG',
            'sequences': [TaPGP, pga, mg]
        },
        'P_TfPGP_L_PGA_C_MG': {
            'name': 'P_TfPGP_L_PGA_C_MG',
            'sequences': [TfPGP, pga, mg]
        },
        'P_CnPGP_L_GOL_C_MG': {
            'name': 'P_CnPGP_L_GOL_C_MG',
            'sequences': [CnPGP, gol, mg]
        },
        'P_TaPGP_L_GOL_C_MG': {
            'name': 'P_TaPGP_L_GOL_C_MG',
            'sequences': [TaPGP, gol, mg]
        },
        'P_TfPGP_L_GOL_C_MG': {
            'name': 'P_TfPGP_L_GOL_C_MG',
            'sequences': [TfPGP, gol, mg]
        }
    }


@pytest.fixture
def one_v_many_ligands_only_results(CnPGP, TaPGP, TfPGP, SePGP, pga, gol):
    return {
        'P_CnPGP_P_SePGP_L_PGA': {
            'name': 'P_CnPGP_L_PGA',
            'sequences': [CnPGP, SePGP, pga]
        },
        'P_TaPGP_P_SePGP_L_PGA': {
            'name': 'P_TaPGP_L_PGA',
            'sequences': [TaPGP, SePGP, pga]
        },
        'P_TfPGP_P_SePGP_L_PGA': {
            'name': 'P_TfPGP_L_PGA',
            'sequences': [TfPGP, SePGP, pga]
        },
        'P_CnPGP_P_SePGP_L_GOL': {
            'name': 'P_CnPGP_L_GOL',
            'sequences': [CnPGP, SePGP, gol]
        },
        'P_TaPGP_P_SePGP_L_GOL': {
            'name': 'P_TaPGP_L_GOL',
            'sequences': [TaPGP, SePGP, gol]
        },
        'P_TfPGP_P_SePGP_L_GOL': {
            'name': 'P_TfPGP_L_GOL',
            'sequences': [TfPGP, SePGP, gol]
        }
    }


@pytest.fixture
def many_v_many_no_ligand_no_cofactor_results(CnPGP, TaPGP, TfPGP):
    return {
        'P_CnPGP_P_TaPGP': {
            'name': 'P_CnPGP_P_TaPGP',
            'sequences': [CnPGP, TaPGP]
        },
        'P_CnPGP_P_TfPGP': {
            'name': 'P_CnPGP_P_TfPGP',
            'sequences': [CnPGP, TfPGP]
        },
        'P_TaPGP_P_TfPGP': {
            'name': 'P_TaPGP_P_TfPGP',
            'sequences': [TaPGP, TfPGP]
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
                              anchor_sequences=anchor_sequences, protein_comparison_type='one_v_many')
    

    assert result == one_v_many_ligands_only_results


def test_generateJSONs_one_v_many_ligands_only(
        fasta_sequences, cofactor_list_empty,
        many_v_many_no_ligand_no_cofactor_results):

    result = af.generateJSONs(fasta_sequences,
                              cofactor_list=cofactor_list_empty, protein_comparison_type='many_v_many')

    assert result == many_v_many_no_ligand_no_cofactor_results
