import ExcelJS from 'exceljs';


export enum MIME_FILE_FORMATS {
    csv = 'text/csv',
    tsv = 'text/tsv',
    xlsx = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
};


export async function exportDataToBlob<T extends Record<string, any>>(
    data: T[],
    fileFormat: MIME_FILE_FORMATS,
): Promise<Blob> {

    if (data.length == 0) {
        throw new Error('cannot export empty data');
    }
    if (!(Object.values(MIME_FILE_FORMATS).includes(fileFormat))) {
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
        case MIME_FILE_FORMATS.csv:
            buffer = await workbook.csv.writeBuffer();
            break;
        case MIME_FILE_FORMATS.tsv:
            buffer = await workbook.csv.writeBuffer({ formatterOptions: { delimiter: '\t' } });
            break;
        case MIME_FILE_FORMATS.xlsx:
            buffer = await workbook.xlsx.writeBuffer();
            break;
    }

    return new Blob([buffer], { type: fileFormat });
}


export enum FILE_EXPORT_STATUS {
    idle = 'Idle',
    inProgress = 'In Progress',
    success = 'Success',
    error = 'Error'
}


export async function handleExportFile<T extends Record<string, any>>(
    data: T[],
    fileFormat: MIME_FILE_FORMATS,
    onChangeStatus: (status: FILE_EXPORT_STATUS) => void,
    fileNamePrefix: string,
) {
    onChangeStatus(FILE_EXPORT_STATUS.inProgress)

    try {
        const blob = await exportDataToBlob(
            data,
            fileFormat,
        );

        const fileFormatEntry = Object.entries(MIME_FILE_FORMATS).find(([_, mimeType]) => fileFormat === mimeType)
        if (fileFormatEntry === undefined) {
            console.error('error exporting file')
            return
        }
        const filename = `${fileNamePrefix}-${Date.now()}.${fileFormatEntry[0]}`;

        const a = document.createElement('a');
        a.setAttribute('href', window.URL.createObjectURL(blob));
        a.setAttribute('download', filename);
        a.click();
        a.remove();

        onChangeStatus(FILE_EXPORT_STATUS.idle)
    } catch (err) {
        console.error(`Error export names file: ${err}`);
        onChangeStatus(FILE_EXPORT_STATUS.error)
    }
};