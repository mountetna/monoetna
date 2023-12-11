import _ from 'lodash';

import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameCompleteCreateName, UNSET_VALUE } from '../../models';
import { NamesState } from './reducer';
import { listToIdObject } from '../../utils/object';
import { RulesStateSliceForCompleteCreateNames } from '../../selectors/global';
import { renderCounter, renderTokens } from '../../utils/names';
import { createLocalId } from '../../utils/models';



export function addCompleteCreateNames(
    completeCreateNames: CompleteCreateName[],
    createNameCompleteCreateNames: CreateNameCompleteCreateName[],
    completeCreateNameParents: CompleteCreateNameParent[],
    state: NamesState,
): NamesState {
    const newCompleteCreateNames = {
        byLocalId: {
            ...state.completeCreateNames.byLocalId,
            ...listToIdObject(completeCreateNames, 'localId'),
        }
    };

    const newCreateNameCompleteCreateNames = {
        byLocalId: { ...state.createNameCompleteCreateNames.byLocalId },
        byCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCreateNameLocalId },
        byCompleteCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCompleteCreateNameLocalId },
    };

    for (const cnccn of createNameCompleteCreateNames) {
        newCreateNameCompleteCreateNames.byLocalId[cnccn.localId] = cnccn;
        newCreateNameCompleteCreateNames.byCreateNameLocalId[cnccn.createNameLocalId] = cnccn.localId;

        if (!(cnccn.completeCreateNameLocalId in newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId)) {
            newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId] = [];
        }
        newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId].push(cnccn.localId);
    }

    const newCompleteCreateNameParents = {
        byLocalId: {
            ...state.completeCreateNameParents.byLocalId,
            ...listToIdObject(completeCreateNameParents, 'localId'),
        },
        byParentLocalId: { ...state.completeCreateNameParents.byParentLocalId },
        byChildLocalId: { ...state.completeCreateNameParents.byChildLocalId },
        byParentLocalIdByChildRenderedTokensByChildCounterValue: _.cloneDeep(
            state.completeCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
        ),
    };

    for (const ccnp of completeCreateNameParents) {
        if (!(ccnp.parentCompleteCreateNameLocalId in newCompleteCreateNameParents.byParentLocalId)) {
            newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId] = [];
        }
        newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId].push(ccnp.localId);

        if (!(ccnp.completeCreateNameLocalId in newCompleteCreateNameParents.byChildLocalId)) {
            newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId] = [];
        }
        newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId].push(ccnp.localId);


        const childCompleteCreateName = state.completeCreateNames.byLocalId[ccnp.completeCreateNameLocalId] = newCompleteCreateNames.byLocalId[ccnp.completeCreateNameLocalId];

        const hierarchyCounterValue = childCompleteCreateName.counterValue != undefined ? childCompleteCreateName.counterValue : UNSET_VALUE;

        if (!(ccnp.parentCompleteCreateNameLocalId in newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue)) {
            newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue[ccnp.parentCompleteCreateNameLocalId] = {};
        }

        const byParentLocalId = newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue[ccnp.parentCompleteCreateNameLocalId];

        if (!(childCompleteCreateName.value in byParentLocalId)) {
            // @ts-ignore
            byParentLocalId[childCompleteCreateName.value] = {};
        }
        byParentLocalId[childCompleteCreateName.value][hierarchyCounterValue] = ccnp.localId;
    }

    return {
        ...state,
        completeCreateNames: newCompleteCreateNames,
        createNameCompleteCreateNames: newCreateNameCompleteCreateNames,
        completeCreateNameParents: newCompleteCreateNameParents,
    };
}

export function removeCompleteCreateNames(
    createNameCompleteCreateNames: CreateNameCompleteCreateName[],
    state: NamesState,
): NamesState {
    // find CompleteCreateNames and Parents from CreateNameCompleteCreateName
    const completeCreateNames: CompleteCreateName[] = [];
    const completeCreateNameParents = new Set<CompleteCreateNameParent>();

    const newCreateNameCompleteCreateNames = {
        byLocalId: { ...state.createNameCompleteCreateNames.byLocalId },
        byCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCreateNameLocalId },
        byCompleteCreateNameLocalId: { ...state.createNameCompleteCreateNames.byCompleteCreateNameLocalId },
    };

    for (const cnccn of createNameCompleteCreateNames) {
        delete newCreateNameCompleteCreateNames.byLocalId[cnccn.localId];
        delete newCreateNameCompleteCreateNames.byCreateNameLocalId[cnccn.createNameLocalId];

        if (cnccn.completeCreateNameLocalId in newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId) {

            newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId] =
                newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId]
                    .filter(candidateCnccnLocalId => candidateCnccnLocalId != cnccn.localId);

            if (newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId].length == 0) {
                delete newCreateNameCompleteCreateNames.byCompleteCreateNameLocalId[cnccn.completeCreateNameLocalId];

                // and delete CompleteCreateName and Parents since there aren't any CreateName associations left
                completeCreateNames.push(state.completeCreateNames.byLocalId[cnccn.completeCreateNameLocalId]);

                const parentsToDelete = [
                    ...(state.completeCreateNameParents.byParentLocalId[cnccn.completeCreateNameLocalId] || []),
                    ...(state.completeCreateNameParents.byChildLocalId[cnccn.completeCreateNameLocalId] || []),
                ];
                parentsToDelete.forEach(ccnpLocalId => {
                    completeCreateNameParents.add(state.completeCreateNameParents.byLocalId[ccnpLocalId]);
                });
            }
        }
    }

    const newCompleteCreateNames = {
        byLocalId: { ...state.completeCreateNames.byLocalId },
    };

    for (const ccn of completeCreateNames) {
        delete newCompleteCreateNames.byLocalId[ccn.localId];
    }

    const newCompleteCreateNameParents = {
        byLocalId: { ...state.completeCreateNameParents.byLocalId },
        byParentLocalId: { ...state.completeCreateNameParents.byParentLocalId },
        byChildLocalId: { ...state.completeCreateNameParents.byChildLocalId },
        byParentLocalIdByChildRenderedTokensByChildCounterValue: _.cloneDeep(
            state.completeCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
        ),
    };

    for (const ccnp of completeCreateNameParents) {
        delete newCompleteCreateNameParents.byLocalId[ccnp.localId];

        if (ccnp.parentCompleteCreateNameLocalId in newCompleteCreateNameParents.byParentLocalId) {

            newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId] =
                newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId]
                    .filter(candidateCcnpLocalId => candidateCcnpLocalId != ccnp.localId);

            if (newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId].length == 0) {
                delete newCompleteCreateNameParents.byParentLocalId[ccnp.parentCompleteCreateNameLocalId];
            }
        }

        if (ccnp.completeCreateNameLocalId in newCompleteCreateNameParents.byChildLocalId) {

            newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId] =
                newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId]
                    .filter(candidateCcnpLocalId => candidateCcnpLocalId != ccnp.localId);

            if (newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId].length == 0) {
                delete newCompleteCreateNameParents.byChildLocalId[ccnp.completeCreateNameLocalId];
            }
        }

        // remove from byParentLocalIdByChildRenderedTokensByChildCounterValue
        const childCompleteCreateName = state.completeCreateNames.byLocalId[ccnp.completeCreateNameLocalId];

        delete newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
        [ccnp.parentCompleteCreateNameLocalId][childCompleteCreateName.value]
        [childCompleteCreateName.counterValue != undefined ? childCompleteCreateName.counterValue : UNSET_VALUE];

        if (
            Object.keys(newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId][childCompleteCreateName.value]).length == 0
        ) {
            delete newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId][childCompleteCreateName.value];
        }

        if (
            Object.keys(newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId]).length == 0
        ) {
            delete newCompleteCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue
            [ccnp.parentCompleteCreateNameLocalId];
        }
    }

    return {
        ...state,
        completeCreateNames: newCompleteCreateNames,
        createNameCompleteCreateNames: newCreateNameCompleteCreateNames,
        completeCreateNameParents: newCompleteCreateNameParents,
    };
}

export function isReplaceOrSearchGroup(createNameGroupLocalId: string, state: NamesState): boolean {
    return (
        state.createNameGroups.replaceLocalIds.has(createNameGroupLocalId)
        || state.createNameGroups.searchLocalIds.has(createNameGroupLocalId)
    );
}

// TODO: handle multiparent rules
// TODO: ID should be monotonically increasing
export function addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
    createNameGroupsLocalIds: string[],
    rulesStateSlice: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    let newState = { ...state };

    // for each CreateNameGroup
    for (const cngLocalId of createNameGroupsLocalIds) {
        const completeCreateNamesToAdd: CompleteCreateName[] = [];
        const createNameCompleteCreateNamesToAdd: CreateNameCompleteCreateName[] = [];
        const createNameCompleteCreateNamesToRemove: CreateNameCompleteCreateName[] = [];
        const completeCreateNameParentsToAdd: CompleteCreateNameParent[] = [];

        if (isReplaceOrSearchGroup(cngLocalId, newState)) {
            continue;
        }

        const createNameGroup = newState.createNameGroups.byLocalId[cngLocalId];
        const primaryRuleName = newState.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId].ruleName;
        const ruleNamesHierarchicalList = rulesStateSlice.ruleNamesHierarchicalListByPrimaryRuleName[primaryRuleName];

        const createNames = newState.createNames.byCreateNameGroupLocalId[cngLocalId]
            .map(cnLocalId => newState.createNames.byLocalId[cnLocalId]);

        const sortedCreateNames: CreateName[] = _.sortBy(createNames, [(cn: CreateName) => {
            return ruleNamesHierarchicalList.indexOf(cn.ruleName);
        }]);

        // keep track of hierarchy completeness
        const completeCreateNames: (CompleteCreateName | undefined)[] = [];
        let completeHierarchy = true;

        // for each CreateName
        for (const [idx, createName] of sortedCreateNames.entries()) {
            const createNameTokenValues = (newState.createNameTokenValues.byCreateNameLocalId[createName.localId] || [])
                .map(cntvLocalId => newState.createNameTokenValues.byLocalId[cntvLocalId]);

            const actualTokenCount = createNameTokenValues.length;
            const targetTokenCount = rulesStateSlice.ruleTokenLocalIdsByRuleName[createName.ruleName].length;
            const counterRequired = rulesStateSlice.counterRequiredByRuleName[createName.ruleName];
            const hasCounterValue = createName.ruleCounterValue != undefined;

            let parentCompleteCreateName: CompleteCreateName | undefined = undefined;

            if (idx > 0) {
                parentCompleteCreateName = completeCreateNames[idx - 1];

                if (parentCompleteCreateName == undefined) {
                    completeHierarchy = false;
                }
            }

            // detect completeness
            if (
                targetTokenCount === actualTokenCount
                && (
                    !counterRequired
                    || hasCounterValue
                )
                && completeHierarchy
            ) {
                // is complete (at this level)

                let completeCreateNameParentLocalId: string | undefined = undefined;
                const hierarchyRenderedTokenValue = renderTokens(
                    createNameTokenValues,
                    rulesStateSlice.ruleTokensByLocalId,
                    rulesStateSlice.tokenValuesByLocalId,
                );
                const hierarchyRenderedCounterValue = createName.ruleCounterValue != undefined ? Number.parseInt(renderCounter(createName)) : UNSET_VALUE;

                // find completeCreateNameParentLocalId
                if (parentCompleteCreateName != undefined || idx == 0) {
                    const parentCompleteCreateNameLocalId = parentCompleteCreateName != undefined ? parentCompleteCreateName.localId : UNSET_VALUE;

                    const renderedTokenValuesDict = newState.completeCreateNameParents
                        .byParentLocalIdByChildRenderedTokensByChildCounterValue[parentCompleteCreateNameLocalId];

                    if (renderedTokenValuesDict != undefined) {
                        const renderedCounterValuesDict = renderedTokenValuesDict[hierarchyRenderedTokenValue];

                        if (renderedCounterValuesDict != undefined) {
                            completeCreateNameParentLocalId = renderedCounterValuesDict[hierarchyRenderedCounterValue];
                        }
                    }
                }

                let completeCreateNameParent: CompleteCreateNameParent | undefined = undefined;

                if (completeCreateNameParentLocalId) {
                    completeCreateNameParent = newState.completeCreateNameParents.byLocalId[completeCreateNameParentLocalId];
                }


                // create CompleteCreateName and CreateName association if doesn't exist
                let completeCreateName: CompleteCreateName | undefined = undefined;

                if (completeCreateNameParent == undefined) {
                    completeCreateName = {
                        localId: createLocalId(),
                        value: renderTokens(
                            createNameTokenValues,
                            rulesStateSlice.ruleTokensByLocalId,
                            rulesStateSlice.tokenValuesByLocalId,
                        ),
                        counterValue: createName.ruleCounterValue,
                    };

                    completeCreateNamesToAdd.push(completeCreateName);

                    // remove existing association and its parents if exist
                    const existingAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId];
                    const existingAssociation = newState.createNameCompleteCreateNames.byLocalId[existingAssociationLocalId];

                    if (existingAssociation) {
                        createNameCompleteCreateNamesToRemove.push(existingAssociation);
                    }

                    // create new association
                    const completeCreateNameAssociation: CreateNameCompleteCreateName = {
                        localId: createLocalId(),
                        createNameLocalId: createName.localId,
                        completeCreateNameLocalId: completeCreateName.localId,
                    };
                    createNameCompleteCreateNamesToAdd.push(completeCreateNameAssociation);
                } else {
                    completeCreateName = newState.completeCreateNames.byLocalId[completeCreateNameParent.completeCreateNameLocalId];
                    const existingAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId];
                    const existingAssociation = newState.createNameCompleteCreateNames.byLocalId[existingAssociationLocalId];

                    if (existingAssociation?.completeCreateNameLocalId != completeCreateName.localId) {
                        // create association
                        const completeCreateNameAssociation: CreateNameCompleteCreateName = {
                            localId: createLocalId(),
                            createNameLocalId: createName.localId,
                            completeCreateNameLocalId: completeCreateName.localId,
                        };
                        createNameCompleteCreateNamesToAdd.push(completeCreateNameAssociation);

                        if (existingAssociation) {
                            createNameCompleteCreateNamesToRemove.push(existingAssociation);
                        }
                    }
                }

                // create CompleteCreateNameParent if doesn't exist
                if (completeCreateNameParent == undefined) {
                    const parent: CompleteCreateNameParent = {
                        localId: createLocalId(),
                        completeCreateNameLocalId: completeCreateName.localId,
                        parentCompleteCreateNameLocalId: parentCompleteCreateName?.localId || UNSET_VALUE,
                    };

                    completeCreateNameParentsToAdd.push(parent);
                }

                completeCreateNames.push(completeCreateName);

            } else {
                // is not complete (at this level)
                // remove all records if exist

                const completeCreateNameAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId];

                if (completeCreateNameAssociationLocalId != undefined) {
                    const completeCreateNameAssociation = newState.createNameCompleteCreateNames.byLocalId[completeCreateNameAssociationLocalId];
                    createNameCompleteCreateNamesToRemove.push(completeCreateNameAssociation);
                }

                completeCreateNames.push(undefined);
            }
        }

        newState = addCompleteCreateNames(
            completeCreateNamesToAdd,
            createNameCompleteCreateNamesToAdd,
            completeCreateNameParentsToAdd,
            removeCompleteCreateNames(createNameCompleteCreateNamesToRemove, newState),
        );
    }

    return newState;
}

// TODO: handle multiparent rules
export function removeCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
    createNameGroupsLocalIds: string[],
    state: NamesState,
): NamesState {
    let newState = { ...state };

    // for each CreateNameGroup
    for (const cngLocalId of createNameGroupsLocalIds) {
        const createNameCompleteCreateNamesToRemove: CreateNameCompleteCreateName[] = [];

        if (isReplaceOrSearchGroup(cngLocalId, newState)) {
            continue;
        }

        const createNames = newState.createNames.byCreateNameGroupLocalId[cngLocalId]
            .map(cnLocalId => newState.createNames.byLocalId[cnLocalId]);

        for (const createName of createNames) {

            const completeCreateNameAssociationLocalId = newState.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId];
            const completeCreateNameAssociation = newState.createNameCompleteCreateNames.byLocalId[completeCreateNameAssociationLocalId];

            if (completeCreateNameAssociation) {
                createNameCompleteCreateNamesToRemove.push(completeCreateNameAssociation);
            }
        }

        newState = removeCompleteCreateNames(createNameCompleteCreateNamesToRemove, newState);
    }

    return newState;
}