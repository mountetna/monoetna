import {stubUrl, cleanStubs, joinedDeferredPromises} from '../../../spec/helpers';
import setupUploadWorker, {
  AddUploadCommand,
  CancelUploadCommand,
  PauseUploadCommand,
  UploadErrorEvent,
  hashBlob,
  Upload,
  ZERO_HASH, UPLOAD_COMPLETE, UPLOAD_ERROR
} from '../uploader';
import {fileKey} from "../../../utils/file";
import {Subscription} from "../../../utils/subscription";

const uploadSubscription = new Subscription();

class TestUploadWorker {
  constructor(opts) {
    const postMessage = this._postMessage = jest.fn();
    const addEventListener = jest.fn();

    this.uploader = setupUploadWorker({postMessage, addEventListener}, opts);
    uploadSubscription.addCleanup(this.uploader.subscription);

    expect(addEventListener).toHaveBeenCalled();
    expect(addEventListener.mock.calls[0][0]).toEqual('message');

    this._messageHandler = addEventListener.mock.calls[0][1];
  }

  dispatch(data) {
    this._messageHandler({data});
  }

  popEvent(predicate) {
    let next;
    do {
      next = this._postMessage.mock.calls.shift();
    } while (next && !predicate(next = next[0]) && (this.prev = next))

    const {prev} = this;
    this.prev = next;
    return [prev, next];
  }
}

function serverUploadResponse({project_name, file_name, author, current_byte_position, next_blob_size, next_blob_hash}) {
  return {
    project_name,
    file_name,
    author,
    current_byte_position,
    next_blob_hash,
    next_blob_size
  }
}

describe('Uploader', () => {
  afterEach(() => {
    uploadSubscription.end();
  });

  const uploadFor = (file) => Upload({
    file,
    file_name: 'test/' + file.name,
    project_name: 'test',
    url: 'http://localhost/' + file.name
  });

  const opts = {maxBlobSize: 10, maxUploads: 2, minBlobSize: 10, debouncerOptions: {maxGating: 0, eager: true}};
  const emptyFile = new File([], 'empty.txt');
  const aFile = new File(['some text would go here'], 'a.txt'); // 24
  const bFile = new File(['different text would go here'], 'b.txt'); // 29
  const shortFile = new File(['hi'], 'short.txt'); // 2

  function stubUploadStart(upload, reset, response = null, status = 200) {
    const {file} = upload;

    let hash;
    if (file.size === 0) {
      hash = Promise.resolve(ZERO_HASH)
    } else {
      hash = hashBlob(file.slice(0, opts.minBlobSize))
    }

    if (response == null) {
      response = serverUploadResponse({
        ...upload,
        next_blob_size: Math.min(file.size, opts.minBlobSize),
      })
    }

    return hash.then(next_blob_hash => {
      return () => stubUrl({
        verb: 'post',
        path: '/' + file.name,
        request: {
          file_size: file.size,
          next_blob_size: Math.min(file.size, opts.minBlobSize),
          next_blob_hash,
          action: 'start',
          reset,
        },
        response,
        status,
        host: 'http://localhost'
      });
    })
  }

  function stubUploadBlob(upload, newBlobSize, response = null, status = 200) {
    const {current_byte_position, file, next_blob_size} = upload;
    const blob_data = file.slice(current_byte_position, current_byte_position + next_blob_size);

    return hashBlob(file.slice(current_byte_position + next_blob_size, current_byte_position + next_blob_size + newBlobSize)).then(newBlobHash => {
      if (response == null) {
        response = serverUploadResponse({
          ...upload,
          next_blob_size: newBlobSize,
          next_blob_hash: newBlobHash,
          current_byte_position: current_byte_position + next_blob_size
        })
      }

      return () => stubUrl({
        verb: 'post',
        path: '/' + file.name,
        request: new FormData({
          action: 'blob',
          blob_data,
          next_blob_size: newBlobSize,
          next_blob_hash: newBlobHash,
          current_byte_position
        }).toString(),
        response,
        status,
        host: 'http://localhost'
      });
    })
  }

  function stubUploadCancel(upload, status = 200, response = {}) {
    const {file, project_name, file_name} = upload;
    return Promise.resolve(() => stubUrl({
      verb: 'post',
      path: '/' + file.name,
      request: {
        project_name,
        file_name,
        action: 'cancel',
      },
      response,
      status,
      host: 'http://localhost'
    }));
  }

  describe('error handling', () => {
    describe('exceptions during blob', () => {
      describe('retryable errors', () => {
        it('retries 500s, 422s, body parsing, and network exceptions up to maxRetries via start, and resets only on 422s', async () => {
          let worker = new TestUploadWorker({...opts, maxUploadFailures: 5, backoffFactor: 1});
          let aUpload = uploadFor(aFile);

          let expectedRequests = await joinedDeferredPromises(
            [
              stubUploadStart(aUpload, false),
              stubUploadBlob({...aUpload, next_blob_size: 10}, 10, {error: 'Oh no'}, 500),
              stubUploadStart({...aUpload, next_blob_size: 10}, false),
              stubUploadBlob({...aUpload, next_blob_size: 10}, 10, (uri, body, cb) => {
                cb(new Error('Network ERR'))
              }),
              stubUploadStart({...aUpload, next_blob_size: 10}, false),
              stubUploadBlob({...aUpload, next_blob_size: 10}, 10, 'Not json!', 200),
              stubUploadStart({...aUpload, next_blob_size: 10}, false),
              stubUploadBlob({...aUpload, next_blob_size: 10}, 10, '', 422),

              // Reset from the above 422
              stubUploadStart(aUpload, true),
              // The final request that will fail due to max retries
              stubUploadBlob({...aUpload, next_blob_size: 10}, 10, {error: 'Last'}, 500),
            ],
          );

          worker.dispatch(AddUploadCommand(aUpload));

          await Promise.all(expectedRequests);
          expect(await worker.uploader.schedule.allPending().catch(e => e.toString())).toEqual('Error: Last')


          let [prev, next] = worker.popEvent(({type}) => type === UPLOAD_ERROR);
          expect([next]).toContainEqual(UploadErrorEvent('error', 'Upload of file test/a.txt', 'Last', {
            ...aUpload,
            next_blob_size: 10,
            status: 'paused'
          }));

          [prev, next] = worker.popEvent(({type}) => type === UPLOAD_ERROR);
          expect(next).toBeUndefined();
        })
      })
    });

    describe('exceptions during start', () => {
      it('pauses the upload and emits an error', async () => {
        let worker = new TestUploadWorker(opts);
        let aUpload = uploadFor(aFile);

        let expectedRequests = await joinedDeferredPromises(
          [stubUploadStart(aUpload, false, {error: 'Error message'}, 500),],
        );

        worker.dispatch(AddUploadCommand(aUpload));

        await worker.uploader.schedule.allPending().catch(e => e);
        await Promise.all(expectedRequests);

        let [prev, next] = worker.popEvent(({uploads} = {}) => uploads[fileKey(aUpload)].status === 'paused');
        expect(next.uploads[fileKey(aUpload)].status).toEqual('paused');

        [prev, next] = worker.popEvent(({type}) => type === UPLOAD_ERROR);
        expect([next]).toContainEqual(UploadErrorEvent('error', 'Upload of file test/a.txt', 'Error message', {
          ...aUpload,
          status: 'paused'
        }));

        [prev, next] = worker.popEvent(({type}) => type === UPLOAD_ERROR);
        expect(next).toBeUndefined();
      })
    })
  });

  describe('prepareNextBlobParams', () => {
    it('passes through a metis_uid if provided', (done) => {
      let worker = new TestUploadWorker(opts);
      let upload = Upload({
        file: aFile,
        file_name: 'test/' + aFile.name,
        project_name: 'test',
        url: 'http://localhost/' + aFile.name,
        metis_uid: '123'
      });

      worker.uploader.prepareNextBlobParams(upload)
      .then((params) => {
        expect(params['metis_uid']).toEqual('123');
        done();
      });
    })

    it('does not include metis_uid if not provided', (done) => {
      let worker = new TestUploadWorker(opts);
      let upload = Upload({
        file: aFile,
        file_name: 'test/' + aFile.name,
        project_name: 'test',
        url: 'http://localhost/' + aFile.name
      });

      worker.uploader.prepareNextBlobParams(upload)
      .then((params) => {
        expect(params.hasOwnProperty('metis_uid')).toEqual(false);
        done();
      });
    })
  });

  describe('queueing and uploading', () => {
    describe('in the face of removing an upload', () => {
      it('stops the cancelled upload, then allows others to proceed', async () => {
        let worker = new TestUploadWorker({...opts, maxUploads: 1});

        let bUpload = uploadFor(bFile);
        let aUpload = uploadFor(aFile);

        // We'll start the b upload, and pause it after getting one chunk up.
        let expectedRequests = await joinedDeferredPromises(
          [
            stubUploadStart(bUpload, false),
            stubUploadBlob({...bUpload, next_blob_size: 10}, 10),
          ],
        );

        worker.dispatch(AddUploadCommand(bUpload));
        worker.dispatch(AddUploadCommand(aUpload));

        await Promise.all(expectedRequests);
        // Empty all events
        worker.popEvent(() => false);

        expectedRequests = await joinedDeferredPromises(
          [
            stubUploadStart(aUpload, false),
            stubUploadBlob({...aUpload, next_blob_size: 10}, 10),
            stubUploadBlob({...aUpload, next_blob_size: 10, current_byte_position: 10}, 3),
            stubUploadBlob({...aUpload, next_blob_size: 3, current_byte_position: 20}, 0),
          ],
          [
            stubUploadCancel(bUpload),
          ]
        );

        worker.dispatch(CancelUploadCommand(bUpload));

        await worker.uploader.schedule.allPending();
        await Promise.all(expectedRequests);

        let [_, next] = worker.popEvent(({uploads}) => uploads);
        expect(next.uploads).not.toHaveProperty(fileKey(bUpload));


        // The removed upload is not 'completed'.
        [_, next] = worker.popEvent(({type}) => type === UPLOAD_COMPLETE);
        expect(fileKey(next.upload)).not.toEqual(fileKey(bUpload));
        [_, next] = worker.popEvent(({type}) => type === UPLOAD_COMPLETE);
        expect(next).toBeUndefined();
      });
    });

    describe('in the face of pausing', () => {
      it('stops, then resumes on add, without yielding to other uploads', async () => {
        let worker = new TestUploadWorker({...opts, maxUploads: 1});
        let bUpload = uploadFor(bFile);
        let aUpload = uploadFor(aFile);

        // We'll start the b upload, and pause it after getting one chunk up.
        let expectedRequests = await joinedDeferredPromises(
          [
            stubUploadStart(bUpload, false),
            stubUploadBlob({...bUpload, next_blob_size: 10}, 10),
          ],
        );

        worker.dispatch(AddUploadCommand(bUpload));
        worker.dispatch(AddUploadCommand(aUpload));

        await Promise.all(expectedRequests);

        worker.dispatch(PauseUploadCommand(bUpload));

        let [_, next] = worker.popEvent(({uploads = {}}) => uploads[fileKey(bUpload)].status === 'paused');
        expect([next.uploads[fileKey(aUpload)].status === 'queued']);


        expectedRequests = await joinedDeferredPromises(
          [
            stubUploadStart(bUpload, false),
            stubUploadBlob({...bUpload, next_blob_size: 10}, 10),
            stubUploadBlob({...bUpload, next_blob_size: 10, current_byte_position: 10}, 8),
            stubUploadBlob({...bUpload, next_blob_size: 8, current_byte_position: 20}, 0),
          ],
          [
            stubUploadStart(aUpload, false),
          ],
        );

        worker.dispatch(AddUploadCommand(bUpload));

        await Promise.all(expectedRequests);

        [_, next] = worker.popEvent(({type}) => type === UPLOAD_COMPLETE);
        expect([fileKey(next.upload)]).toContainEqual(fileKey(bUpload));
      })
    });

    it('allows up to maxUploads, and starts Uploads as others finish', async () => {
      let worker = new TestUploadWorker(opts);

      let emptyUpload = uploadFor(emptyFile);
      let aUpload = uploadFor(aFile);
      let bUpload = uploadFor(bFile);
      let shortUpload = uploadFor(shortFile);

      let expectedRequests = await joinedDeferredPromises(
        [stubUploadStart(emptyUpload, false), stubUploadBlob(emptyUpload, 0)],
        [
          stubUploadStart(bUpload, false),
          stubUploadBlob({...bUpload, next_blob_size: 10}, 10),
          stubUploadBlob({...bUpload, next_blob_size: 10, current_byte_position: 10}, 8),
          stubUploadBlob({...bUpload, next_blob_size: 8, current_byte_position: 20}, 0),
        ],
        [
          stubUploadStart(aUpload, false),
          stubUploadBlob({...aUpload, next_blob_size: 10}, 10),
          stubUploadBlob({...aUpload, next_blob_size: 10, current_byte_position: 10}, 3),
          stubUploadBlob({...aUpload, next_blob_size: 3, current_byte_position: 20}, 0),
        ],
        [
          stubUploadStart(shortUpload, false),
          stubUploadBlob({...shortUpload, next_blob_size: 2}, 0),
        ],
      );

      worker.dispatch(AddUploadCommand(emptyUpload));
      worker.dispatch(AddUploadCommand(aUpload));
      worker.dispatch(AddUploadCommand(bUpload));
      worker.dispatch(AddUploadCommand(shortUpload));

      await Promise.all(expectedRequests);
      await worker.uploader.schedule.allPending();

      let [prev, next] = worker.popEvent(({uploads}) => Object.keys(uploads || {}).length === 4);

      expect([next.uploads[fileKey(emptyUpload)].status]).toContainEqual('active');
      expect([next.uploads[fileKey(aUpload)].status]).toContainEqual('active');
      expect([next.uploads[fileKey(bUpload)].status]).toContainEqual('queued');
      expect([next.uploads[fileKey(shortUpload)].status]).toContainEqual('queued');

      [prev, next] = worker.popEvent(({type}) => type === UPLOAD_COMPLETE);
      expect([prev.uploads[fileKey(emptyUpload)].status]).toContainEqual('complete');
      expect([fileKey(next.upload)]).toContainEqual(fileKey(emptyUpload));

      [prev, next] = worker.popEvent(({type}) => type === UPLOAD_COMPLETE);
      expect([prev.uploads[fileKey(aUpload)].status]).toContainEqual('complete');

      // Now that two have completed, the others should be active and running
      expect([prev.uploads[fileKey(bUpload)].status]).toContainEqual('active');
      expect([prev.uploads[fileKey(shortUpload)].status]).toContainEqual('active');
      expect([fileKey(next.upload)]).toContainEqual(fileKey(aUpload));

      [prev, next] = worker.popEvent(({type}) => type === UPLOAD_COMPLETE);
      expect([prev.uploads[fileKey(bUpload)].status]).toContainEqual('complete');
      expect([fileKey(next.upload)]).toContainEqual(fileKey(bUpload));

      // Verify upload_speeds are being processed
      expect([prev.uploads[fileKey(bUpload)].upload_speeds.length]).toContainEqual(3);
    });
  });
});
