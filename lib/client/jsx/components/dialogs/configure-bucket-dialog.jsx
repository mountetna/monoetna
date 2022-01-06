import * as React from 'react';
import { connect } from 'react-redux';

import ConfigRow from './config-row';

const VALID_EMAIL = /^(([^<>()\[\]\\.,;:\s@"]+(\.[^<>()\[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;

const VALID_NAME = /^\w+$/;

class ConfigureBucketDialog extends React.Component {
  constructor(props) {
    super(props)
    this.state = {
      bucket_name: '',
      description: '',
      access: '',
      access_email: '',
      errors: []
    }
  }

  componentDidMount() {
    let { bucket_name='', description='', access='' } = this.props;
    let access_email = '';

    if (/\@/.test(access)) {
      access_email = access;
      access = '';
    }

    this.setState({ bucket_name: bucket_name || '',
      description: description || '',
      access: access || '',
      access_email: access_email || '' });
  }

  updateValue(attribute, e) {
    this.setState({ [attribute]: e.target.value == 'none' ? null : e.target.value });
  }

  submit() {
    let { createBucket, updateBucket, dismissDialog } = this.props;
    let { bucket_name, description, access, access_email } = this.state;

    let errors = [];

    if (!VALID_NAME.test(bucket_name)) errors.push('Invalid bucket name.')

    if (access_email) {
      access = access_email.split(/[\s\;\,]+/)
      if (access.some(email => !VALID_EMAIL.test(email))) {
        errors.push('Invalid email ids in access list.');
      }
      else access = access.join(',');
    }

    if (!access) errors.push('Access level must be set');

    if (errors.length > 0) {
      this.setState({ errors });
      return;
    }

    if (createBucket)
      createBucket({ bucket_name, description, access, owner: 'metis' });
    else {
      updateBucket({
        bucket_name: this.props.bucket_name,
        new_bucket_name: bucket_name == this.props.bucket_name ? undefined : bucket_name,
        description, access
      });
    }
    dismissDialog();
  }

  render() {
    let { bucket_name, description, access, access_email, errors } = this.state;
    let { createBucket, onClick } = this.props;
    return <div className='config-dialog'>
      <div className='title'>Bucket Configuration</div>
      <div className='errors'>{ errors.map(error => <span className='error'>{error}</span>) }</div>
      <ConfigRow label='Name'>
        <input type='text' pattern="\w+" placeholder='E.g. bucket_name' value={ bucket_name } onChange={ this.updateValue.bind(this, 'bucket_name') } />
      </ConfigRow>
      <ConfigRow label='Description'>
        <input type='text' value={ description } onChange={ this.updateValue.bind(this, 'description') } />
      </ConfigRow>
      <ConfigRow label='Access'>
        {
          !access_email && <select value={ access || 'none' } onChange={ this.updateValue.bind(this, 'access') }>
          {
            [ null, 'administrator', 'editor', 'viewer' ].map( v =>
              <option key={v || 'none'} value={v || 'none'}>
                {v || 'Select level'}
              </option>
            )
          }
          </select>
        }
        { !access_email && !access && <div>or</div> }
        {
          !access && <input disabled={ access } type='text' placeholder='List email ids' value={ access_email } onChange={ this.updateValue.bind(this, 'access_email') } />
        }
      </ConfigRow>
      <div className='submit'>
        <span className='button' disabled={ !(access || access_email) || !bucket_name } onClick={ this.submit.bind(this) }>
          { createBucket ? 'Create' : 'Update' }
        </span>
      </div>
    </div>;
  }
}

export default connect(
  null,
  (dispatch) => ({
    dismissDialog: () => dispatch({type:'DISMISS_DIALOG'})
  })
)(ConfigureBucketDialog);
