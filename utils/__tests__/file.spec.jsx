import { fileKey, folderKey, filePath } from '../file';

describe('File Util', () => {
  it('generates a file key from an object', () => {
    const result = fileKey({
      project_name: 'labors',
      file_name: 'lion.txt'
    });

    expect(result).toEqual('labors:lion.txt');
  });

  it('generates a folder key from an object', () => {
    const result = folderKey({
      project_name: 'labors',
      folder_name: 'Nemean Lion/weaknesses'
    });

    expect(result).toEqual('labors:Nemean Lion/weaknesses');
  });

  it('generates a file path with null folder_name', () => {
    const result = filePath(null, 'stats.txt');

    expect(result).toEqual('stats.txt');
  });

  it('generates a file path with null file_name', () => {
    const result = filePath('blueprints', null);

    expect(result).toEqual('blueprints/');
  });

  it('generates a file path with null file_name and null folder_name', () => {
    const result = filePath(null, null);

    expect(result).toEqual('');
  });

  it('generates a file path given file_name and folder_name', () => {
    const result = filePath('blueprints/helmet', 'helmet.jpg');

    expect(result).toEqual('blueprints/helmet/helmet.jpg');
  });
});
