def pandas_to_latex_table(df, row_name, col_name, file_path):
    col_headers = list(df.columns)
    row_headers = list(df.index)
    col_size = len(col_headers)

    table_str_full = ''

    table_str_begin = '\\begin{table}[H]\n' \
                      '\\centering' \
                      '\\captionsetup{width=.75\\textwidth}'
    table_str_full += table_str_begin + '\n'

    table_str_tabular = '\\begin{tabular}{|p{2cm}|'
    for i in range(col_size):
        table_str_tabular += 'p{2cm}|'
    table_str_tabular += '}'
    table_str_full += f'{table_str_tabular}\n\\hline\n'

    table_str_header = '\\textbf{${0}$/${1}$}'.format(str(row_name), str(col_name))
    for _ in col_headers:
        table_str_header += f' & \\textbf{{0}}'
    table_str_header += ' \\\\ \\hline'
    table_str_full += table_str_header + '\n'

    for row_header in row_headers:
        table_str_row_value = f'\\textbf{{0}}'
        for col_header in col_headers:
            table_str_row_value += f' & {str(df[col_header][row_header])}'
        table_str_row_value += ' \\\\ \\hline \n'
        table_str_full += table_str_row_value + '\n'

    table_str_footer = '\\end{tabular}\n' \
                       '\\caption{}\n' \
                       '\\label{tab:}\n' \
                       '\\end{table}'
    table_str_full += table_str_footer + '\n'

    with open(file_path, 'w') as out_file:
        out_file.write(table_str_full)
    return table_str_full
