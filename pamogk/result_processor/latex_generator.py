def pandas_to_latex_table(df, row_name, col_name, file_path):
    col_headers = list(df.columns)
    row_headers = list(df.index)
    col_size = len(col_headers)

    # begin
    table_str_full = '\\begin{table}[H]\n\\centering\\captionsetup{width=.75\\textwidth}\n'

    # col setup
    table_str_full += '\\begin{tabular}{|p{2cm}|' + 'p{2cm}|' * col_size + '}\n\\hline\n'
    table_str_full += f'\\textbf{{${row_name}$/${col_name}$}}' + ' & \\textbf{0}' * len(col_headers) + ' \\\\ \\hline\n'

    # rows
    for row_header in row_headers:
        table_str_row_value = ' & '.join(str(df[col_header][row_header]) for col_header in col_headers)
        table_str_full += '\\textbf{0} & ' + table_str_row_value + ' \\\\ \\hline \n'

    # end
    table_str_full += '\\end{tabular}\n\\caption{}\n\\label{tab:}\n\\end{table}\n'

    with open(file_path, 'w') as out_file:
        out_file.write(table_str_full)
    return table_str_full
