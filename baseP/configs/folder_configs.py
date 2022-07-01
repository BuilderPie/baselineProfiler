import os
def formOutput(output):
    for folder in [
       #figure and table structrue for CRISPR
       os.path.join(output, 'CCLE_evaluator'),
       os.path.join(output, 'CCLE_evaluator', 'tables'),
       os.path.join(output, 'CCLE_evaluator', 'tables', 'JSON'),
       os.path.join(output, 'CCLE_evaluator', 'Exprsn'),
       os.path.join(output, 'CCLE_evaluator', 'Exprsn', 'tables'),
       os.path.join(output, 'CCLE_evaluator', 'Exprsn', 'plots'),
       os.path.join(output, 'CCLE_evaluator', 'Proteomics'),
       os.path.join(output, 'CCLE_evaluator', 'Proteomics', 'tables'),
       os.path.join(output, 'CCLE_evaluator', 'Proteomics', 'plots'),
       os.path.join(output, 'CCLE_evaluator', 'CRISPR_Broad'),
       os.path.join(output, 'CCLE_evaluator', 'CRISPR_Broad', 'tables'),
       os.path.join(output, 'CCLE_evaluator', 'CRISPR_Broad', 'plots'),

       os.path.join(output, 'GTEx_evaluator'),
       os.path.join(output, 'GTEx_evaluator', 'Exprsn'),
       os.path.join(output, 'GTEx_evaluator', 'Exprsn', 'tables'),
       os.path.join(output, 'GTEx_evaluator', 'Exprsn', 'plots'),

       os.path.join(output, 'HPA_evaluator'),
       os.path.join(output, 'HPA_evaluator', 'Exprsn'),
       os.path.join(output, 'HPA_evaluator', 'Exprsn', 'tables'),
       os.path.join(output, 'HPA_evaluator', 'Exprsn', 'plots'),
       
       
       ### html report
        os.path.join(output, 'html'),
        os.path.join(output, 'html', 'figs'),
        os.path.join(output, 'html', 'figs', 'CCLE_evaluator'),
        os.path.join(output, 'html', 'figs', 'CCLE_evaluator', 'Exprsn'),
        os.path.join(output, 'html', 'figs', 'CCLE_evaluator', 'Proteomics'),
        os.path.join(output, 'html', 'figs', 'CCLE_evaluator', 'CRISPR_Broad'),

        os.path.join(output, 'html', 'figs', 'GTEx_evaluator'),
        os.path.join(output, 'html', 'figs', 'GTEx_evaluator', 'Exprsn'),

        os.path.join(output, 'meta'),
        os.path.join(output, 'meta', 'tables'),
        os.path.join(output, 'meta', 'tables', 'JSON'),
    ]:
        if not os.path.isdir(folder):
            os.makedirs(folder)

