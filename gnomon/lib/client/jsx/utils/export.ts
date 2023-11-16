import ExcelJS from 'exceljs';



export const FILE_FORMATS_TO_MIME = {
    'csv': 'text/csv',
    'tsv': 'text/tsv',
    'xlsx': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
};

export type FileFormat = keyof typeof FILE_FORMATS_TO_MIME


export async function exportDataToBlob<T extends Record<string, any>>(
    data: T[],
    fileFormat: keyof typeof FILE_FORMATS_TO_MIME
): Promise<Blob> {

    const workbook = new ExcelJS.Workbook();
    const sheet = workbook.addWorksheet();

    if (data.length == 0) {
        throw new Error('cannot export empty data');
    }

    // set headers
    sheet.columns = Object.keys(data[0]).map(key => ({
        header: key,
        key,
    }));

    // add data
    sheet.addRows(data);

    let buffer: ExcelJS.Buffer;

    switch (fileFormat) {
        case 'csv':
            buffer = await workbook.csv.writeBuffer();
        case 'tsv':
            buffer = await workbook.csv.writeBuffer({ formatterOptions: { delimiter: '\t' } });
            break;
        case 'xlsx':
            buffer = await workbook.xlsx.writeBuffer();
            break;
        default:
            throw new Error(`unsupported file format: ${fileFormat}`);
    }

    return new Blob([buffer], { type: FILE_FORMATS_TO_MIME[fileFormat] });
}