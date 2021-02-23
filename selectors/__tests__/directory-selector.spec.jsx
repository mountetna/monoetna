import {
  selectFiles,
  selectFolders,
  selectBuckets,
  selectUpload,
  selectUploads
} from '../directory-selector';

describe('Directory Selector', () => {
  it('returns files', () => {
    const files = [{ name: 'stats.txt' }];
    const result = selectFiles({
      directory: { files }
    });

    expect(result).toEqual(files);
  });

  it('returns folders', () => {
    const folders = [{ name: 'blueprints' }];
    const result = selectFolders({
      directory: { folders }
    });

    expect(result).toEqual(folders);
  });

  it('returns buckets', () => {
    const buckets = [{ name: 'learn' }];
    const result = selectBuckets({
      directory: { buckets }
    });

    expect(result).toEqual(buckets);
  });

  it('returns uploads', () => {
    const uploads = { 'athena:wisdom.txt': 'wisdom.txt' };
    const result = selectUploads({
      directory: { uploads }
    });

    expect(result).toEqual(uploads);
  });

  it('returns a single upload', () => {
    const upload = {
      file_name: 'helmet.jpg',
      project_name: 'athena'
    };
    const uploads = {
      'athena:wisdom.txt': 'wisdom.txt',
      'athena:helmet.jpg': 'helmet.jpg'
    };
    const result = selectUpload(
      {
        directory: { uploads }
      },
      upload
    );

    expect(result).toEqual('helmet.jpg');
  });
});
