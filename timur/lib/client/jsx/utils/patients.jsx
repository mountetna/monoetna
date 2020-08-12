// ADHOC labeling of waiver / consent tables.
// TODO: this kind of meta labeling be provided by server.  For now these are restricted to just specific
// columns and projects.

export function isWaiverPatientModel(modelName) {
  return TIMUR_CONFIG.project_name === 'mvir1' && model_name === 'patient';
}

export function documentOnWaiver(document) {
  return document.consent === 'Initial Waiver';
}

export function isWaiverPatientAttribute(modelName, attribute) {
  return isWaiverPatientModel(modelName) && attribute === 'consent';
}
