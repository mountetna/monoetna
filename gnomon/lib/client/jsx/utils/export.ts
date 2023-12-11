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

    if (data.length == 0) {
        throw new Error('cannot export empty data');
    }
    if (!(fileFormat in FILE_FORMATS_TO_MIME)) {
        throw new Error(`unsupported file format: ${fileFormat}`);
    }

    const workbook = new ExcelJS.Workbook();
    const sheet = workbook.addWorksheet();

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
            break;
        case 'tsv':
            buffer = await workbook.csv.writeBuffer({ formatterOptions: { delimiter: '\t' } });
            break;
        case 'xlsx':
            buffer = await workbook.xlsx.writeBuffer();
            break;
    }

    return new Blob([buffer], { type: FILE_FORMATS_TO_MIME[fileFormat] });
}