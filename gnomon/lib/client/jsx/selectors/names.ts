import { createSelector } from 'reselect';
import _ from 'lodash';

import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameCompleteCreateName, CreateNameGroup, CreateNameTokenValue, UNSET_VALUE } from '../models';
import { CompleteCreateNameParentsByRenderedValues, MagmaCheckDuplicateNameRequestState, NamesCreationRequestState, NamesListRequestState, NamesState, namesReducer } from '../reducers/names/reducer';
import { defaultDict } from '../utils/object';
import { State } from '../store';
import { selectRuleNamesHierarchicalListByPrimaryRuleName } from './rules';
import { SearchReplaceCriteria, createSearchReplaceCriteriaFromGroups } from '../utils/names';
import { union } from '../utils/set';
import { Status } from '../utils/models';



export const selectNamesState = (state: State): NamesState => state.names;

export const selectCreateNameGroupsByLocalId = (state: State): Record<string, CreateNameGroup> => state.names.createNameGroups.byLocalId;

export const selectCreateNameGroupWithLocalId = (localId: string) => (state: State): CreateNameGroup => {
    return state.names.createNameGroups.byLocalId[localId];
};

export const selectCreateNameGroupsWithLocalIds = (localIds: string[]) => (state: State): CreateNameGroup[] => {
    return localIds.map(localId => state.names.createNameGroups.byLocalId[localId]);
};

export const selectCreateNameGroupIdsByPrimaryRule: (state: State, omitSearchAndReplaceIds?: boolean, respectFilter?: boolean) => Record<string, string[]> = createSelector(
    [
        (state: State) => state.names,
        (state: State, omitSearchAndReplaceIds: boolean = true) => omitSearchAndReplaceIds,
        (state: State, omitSearchAndReplaceIds: boolean = true, respectFilter: boolean = true) => respectFilter,
    ],
    (names: NamesState, omitSearchAndReplaceIds: boolean, respectFilter: boolean): Record<string, string[]> => {

        const result = defaultDict<string, string[]>(_ => []);

        Object.values(names.createNameGroups.byLocalId).forEach(createNameGroup => {
            if (
                omitSearchAndReplaceIds
                && (
                    names.createNameGroups.searchLocalIds.has(createNameGroup.localId)
                    || names.createNameGroups.replaceLocalIds.has(createNameGroup.localId)
                )
            ) {
                return;
            }
            // only ignore non-filter IDs if there are any in filter
            if (
                respectFilter
                && names.createNameGroups.filterEnabled === true
                && !names.createNameGroups.filterLocalIds.has(createNameGroup.localId)) {
                return;
            }
            const createName: CreateName = names.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId];
            result[createName.ruleName].push(createNameGroup.localId);
        });

        return result;
    }
);

export const selectFilterStatus = (state: State): boolean => state.names.createNameGroups.filterEnabled;

export const selectCreateNamesByLocalId = (state: State): Record<string, CreateName> => state.names.createNames.byLocalId;

export const selectCreateNameWithLocalId = (localId: string) => (state: State): CreateName => {
    return state.names.createNames.byLocalId[localId];
};

export const selectCreateNamesWithLocalIds = (localIds: string[]) => (state: State): CreateName[] => {
    return localIds.map(localId => state.names.createNames.byLocalId[localId]);
};

export const selectCreateNameLocalIdsByGroupId = (state: State): Record<string, string[]> => state.names.createNames.byCreateNameGroupLocalId;

export const selectCreateNameLocalIdsWithGroupId = (groupId: string) => (state: State): string[] => state.names.createNames.byCreateNameGroupLocalId[groupId];

export const selectSelectedCreateNameGroupIds = (state: State): Set<string> => state.names.createNameGroups.selectionLocalIds;

export const selectFilterCreateNameGroupIds = (state: State): Set<string> => state.names.createNameGroups.filterLocalIds;

export const selectFilterEnabledStatus = (state: State): boolean => state.names.createNameGroups.filterEnabled;

export const selectSearchCreateNameGroupIds = (state: State): Set<string> => state.names.createNameGroups.searchLocalIds;

export const selectReplaceCreateNameGroupIds = (state: State): Set<string> => state.names.createNameGroups.replaceLocalIds;

export const selectCreateNameTokenValuesByLocalId = (state: State): Record<string, CreateNameTokenValue> => state.names.createNameTokenValues.byLocalId;

export const selectCreateNameTokenValueLocalIdsByCreateNameLocalId = (state: State): Record<string, string[]> => {
    return state.names.createNameTokenValues.byCreateNameLocalId;
};

export const selectCreateNameTokenValueLocalIdsWithCreateNameLocalId = (createNameLocalId: string) => (state: State): string[] => {
    return state.names.createNameTokenValues.byCreateNameLocalId[createNameLocalId];
};

export const selectCompleteCreateNameParentLocalIdsByRenderedValues = (state: State): CompleteCreateNameParentsByRenderedValues => {
    return state.names.completeCreateNameParents.byParentLocalIdByChildRenderedTokensByChildCounterValue;
};

export const selectRenderedCompleteCreateNamesByLocalId: (state: State) => Record<string, string> = createSelector(
    [
        (state: State) => state,
    ],
    (state: State): Record<string, string> => {

        const renderedNamesByLocalId: Record<string, string> = {};
        const namesState = state.names;

        for (const ccn of Object.values(namesState.completeCreateNames.byLocalId)) {

            let renderedValues: string[] = [];
            const localIdsToProcess: string[] = [ccn.localId];

            while (localIdsToProcess.length) {
                const completeCnLocalId = localIdsToProcess.pop();
                if (!completeCnLocalId) {
                    break;
                }

                const completeCn = namesState.completeCreateNames.byLocalId[completeCnLocalId];

                renderedValues.unshift(
                    completeCn.value
                    + (completeCn.counterValue != undefined ? String(completeCn.counterValue) : '')
                );

                // process parent complete names
                if (completeCn.localId in namesState.completeCreateNameParents.byChildLocalId) {

                    for (const ccnpLocalId of [...namesState.completeCreateNameParents.byChildLocalId[completeCn.localId]]) {
                        const parentCcnLocalId = namesState.completeCreateNameParents.byLocalId[ccnpLocalId].parentCompleteCreateNameLocalId;

                        if (parentCcnLocalId == UNSET_VALUE) {
                            continue;
                        } else if (parentCcnLocalId in renderedNamesByLocalId) {
                            // use precomputed value if exists
                            renderedValues.unshift(renderedNamesByLocalId[parentCcnLocalId]);
                        } else {
                            // queue to compute
                            localIdsToProcess.unshift(parentCcnLocalId);
                        }
                    }
                }
            }

            const renderedValue = renderedValues.reduce((prev, curr) => prev + curr, '');
            renderedNamesByLocalId[ccn.localId] = renderedValue;
        }

        return renderedNamesByLocalId;
    }
);

export const selectRenderedCompleteCreateNamesByCreateNameGroupLocalId: (state: State) => Record<string, string> = createSelector(
    [
        (state: State) => state
    ],
    (state: State): Record<string, string> => {

        const renderedNamesByGroupLocalId: Record<string, string> = {};
        const renderedCompleteCreateNamesByLocalId = selectRenderedCompleteCreateNamesByLocalId(state);

        for (const createNameGroup of Object.values(state.names.createNameGroups.byLocalId)) {

            const targetCompleteCreateNameLocalId: string | undefined = state.names.createNameCompleteCreateNames
                .byCreateNameLocalId[createNameGroup.primaryCreateNameLocalId];

            if (targetCompleteCreateNameLocalId != undefined) {
                renderedNamesByGroupLocalId[createNameGroup.localId] = renderedCompleteCreateNamesByLocalId[targetCompleteCreateNameLocalId];
            }
        }

        return renderedNamesByGroupLocalId;
    }
);

// TODO: handle multiple parents
export const selectSortedCompleteCreateNamesWithCreateNameGroupLocalId = (createNameGroupLocalId: string) => (
    state: State
): (CompleteCreateName | undefined)[] => {

    const ruleNamesHierarchicalListByPrimaryRuleName = selectRuleNamesHierarchicalListByPrimaryRuleName(state);
    const createNameGroup = state.names.createNameGroups.byLocalId[createNameGroupLocalId];
    const createNames = state.names.createNames.byCreateNameGroupLocalId[createNameGroupLocalId]
        .map(createNameLocalId => state.names.createNames.byLocalId[createNameLocalId]);
    const primaryCreateName = state.names.createNames.byLocalId[createNameGroup.primaryCreateNameLocalId];

    return ruleNamesHierarchicalListByPrimaryRuleName[primaryCreateName.ruleName].map(ruleName => {

        const createName = createNames.find(cn => cn.ruleName == ruleName);
        if (createName == undefined) {
            throw new Error(`can't find CreateName where ruleName == ${ruleName}`);
        }

        const completeCreateNameAssociationLocalId = state.names.createNameCompleteCreateNames.byCreateNameLocalId[createName.localId];
        const completeCreateNameAssociation: CreateNameCompleteCreateName | undefined = state.names.createNameCompleteCreateNames
            .byLocalId[completeCreateNameAssociationLocalId];

        if (completeCreateNameAssociation == undefined) {
            return undefined;
        }
        return state.names.completeCreateNames.byLocalId[completeCreateNameAssociation.completeCreateNameLocalId];
    });
};


export interface CompleteCreateNameRequestPayload {
    localId: string
    renderedName: string
    ruleName: string
    // names that must exist for primary CreateNames to be created
    // (i.e., parents of primary CreateNames)
    implicit: boolean
}

export const selectCreateNameCompleteCreateNameLocalIdsByCreateNameLocalId = (state: State): Record<string, string> => {
    return state.names.createNameCompleteCreateNames.byCreateNameLocalId;
};

export const selectCreateNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId = (state: State): Record<string, string[]> => {
    return state.names.createNameCompleteCreateNames.byCompleteCreateNameLocalId;
};

export const selectCreateNameCompleteCreateNamesByLocalId = (state: State): Record<string, CreateNameCompleteCreateName> => {
    return state.names.createNameCompleteCreateNames.byLocalId;
};

export const selectCompleteCreateNameParentLocalIdsByChildLocalId = (state: State): Record<string, string[]> => {
    return state.names.completeCreateNameParents.byChildLocalId;
};

export const selectCompleteCreateNameParentsByLocalId = (state: State): Record<string, CompleteCreateNameParent> => {
    return state.names.completeCreateNameParents.byLocalId;
};

export const selectPrimaryCreateNameCountByCompleteCreateNameLocalId: (state: State) => Record<string, number> = createSelector(
    [(state: State) => state],
    (state: State): Record<string, number> => {

        const primaryCnCountByCcnLocalId: Record<string, number> = {};

        for (const ccnLocalId of Object.keys(state.names.completeCreateNames.byLocalId)) {

            let count = 0;
            const cnccnLocalIds = (state.names.createNameCompleteCreateNames.byCompleteCreateNameLocalId[ccnLocalId] || []);

            for (const cnccnLocalId of cnccnLocalIds) {
                const createName = state.names.createNames.byLocalId[
                    state.names.createNameCompleteCreateNames.byLocalId[cnccnLocalId].createNameLocalId
                ];

                if (state.names.createNameGroups.byLocalId[createName.createNameGroupLocalId].primaryCreateNameLocalId == createName.localId) {
                    count += 1;
                }
            }

            primaryCnCountByCcnLocalId[ccnLocalId] = count;
        }

        return primaryCnCountByCcnLocalId;
    }
);

export const selectPrimaryCreateNameCountWithCompleteCreateNameLocalId = (state: State, completeCreateNameLocalId: string): number => {
    return selectPrimaryCreateNameCountByCompleteCreateNameLocalId(state)[completeCreateNameLocalId];
};

export const selectNamesCreationRequestState = (state: State): NamesCreationRequestState => {
    return state.names.creationRequest;
};

export const selectCheckDuplicateNameRequestWithCreateNameGroupLocalId = (state: State, createNameGroupLocalId: string): MagmaCheckDuplicateNameRequestState | undefined => {
    return state.names.createNameGroups.magmaCheckDuplicateNameRequestsByLocalId[createNameGroupLocalId];
};

export const selectHasMagmaDuplicateWithCreateNameGroupLocalId = (state: State, createNameGroupLocalId: string): boolean => {
    const magmaDuplicateRequestState = state.names.createNameGroups.magmaCheckDuplicateNameRequestsByLocalId[createNameGroupLocalId];
    return magmaDuplicateRequestState && magmaDuplicateRequestState.hasDuplicate || false;
};

export const selectComposeErrorsByCreateNameGroupLocalId: (state: State) => Record<string, boolean> = createSelector(
    [(state: State) => state],
    (state: State): Record<string, boolean> => {
        const composeErrorsByCngLocalId: Record<string, boolean> = {};

        for (const cngLocalId of Object.keys(state.names.createNameGroups.byLocalId)) {

            const primaryCnLocalId = state.names.createNameGroups.byLocalId[cngLocalId].primaryCreateNameLocalId;
            const cnccnLocalId = state.names.createNameCompleteCreateNames.byCreateNameLocalId[primaryCnLocalId];

            let hasLocalDuplicate = false;
            if (cnccnLocalId) {
                const primaryCnCount = selectPrimaryCreateNameCountWithCompleteCreateNameLocalId(
                    state,
                    state.names.createNameCompleteCreateNames.byLocalId[cnccnLocalId].completeCreateNameLocalId,
                );

                hasLocalDuplicate = primaryCnCount > 1;
            }

            composeErrorsByCngLocalId[cngLocalId] = hasLocalDuplicate || selectHasMagmaDuplicateWithCreateNameGroupLocalId(state, cngLocalId);
        }

        return composeErrorsByCngLocalId;
    }
);

export const selectComposeErrorCount = (state: State): number => {
    return Object.values(selectComposeErrorsByCreateNameGroupLocalId(state)).filter(val => val).length;
};

export const selectMagmaNamesListsByRuleName = (state: State): Record<string, NamesListRequestState> => {
    return state.names.magmaNamesListRequestsByRuleName;
};

export const selectSearchVisible = (state: State): boolean => state.names.createNameGroups.searchVisible;

export const selectReplaceVisible = (state: State): boolean => state.names.createNameGroups.replaceVisible;

export const selectSearchReplaceCriteriaFromSearchGroups: (state: State) => SearchReplaceCriteria[] = createSelector(
    [(state: State) => state],
    (state: State): SearchReplaceCriteria[] => createSearchReplaceCriteriaFromGroups(state.names, state.names.createNameGroups.searchLocalIds)
);

export const selectSearchReplaceCriteriaFromReplaceGroups: (state: State) => SearchReplaceCriteria[] = createSelector(
    [(state: State) => state],
    (state: State): SearchReplaceCriteria[] => createSearchReplaceCriteriaFromGroups(state.names, state.names.createNameGroups.replaceLocalIds)
);

export const selectCompleteCreateNamesCreationPayloads: (state: State) => CompleteCreateNameRequestPayload[] = createSelector(
    [
        (state: State) => state
    ],
    (state: State): CompleteCreateNameRequestPayload[] => {

        const completeCreateNameGroupLocalIds = new Set(Object.keys(selectRenderedCompleteCreateNamesByCreateNameGroupLocalId(state)));
        const primaryCreateNameLocalIds = new Set<string>(
            Object.values(state.names.createNameGroups.byLocalId).map(cng => cng.primaryCreateNameLocalId)
        );
        const renderedCompleteCreateNamesByLocalId = selectRenderedCompleteCreateNamesByLocalId(state);
        const completeCreateNamesCreationPayloads: CompleteCreateNameRequestPayload[] = [];

        for (const [ccnLocalId, rendered] of Object.entries(renderedCompleteCreateNamesByLocalId)) {

            let createName: CreateName | undefined = undefined;

            for (const associationLocalId of state.names.createNameCompleteCreateNames.byCompleteCreateNameLocalId[ccnLocalId]) {
                const _createName = state.names.createNames.byLocalId[
                    state.names.createNameCompleteCreateNames.byLocalId[associationLocalId].createNameLocalId
                ];

                // check if belongs to a completed (and error-less) CreateNameGroup
                if (
                    completeCreateNameGroupLocalIds.has(_createName.createNameGroupLocalId)
                    && !selectComposeErrorsByCreateNameGroupLocalId(state)[_createName.createNameGroupLocalId]
                ) {
                    createName = _createName;

                    // break if it's a primary CreateName,
                    // else keep the implicit CreateName association while we keep looking for a primary
                    if (primaryCreateNameLocalIds.has(createName.localId)) {
                        break;
                    }
                }
            }

            if (!createName) {
                continue;
            }

            completeCreateNamesCreationPayloads.push({
                localId: ccnLocalId,
                renderedName: rendered,
                ruleName: createName.ruleName,
                // if it's not a primary CreateName, it's implicit
                implicit: !primaryCreateNameLocalIds.has(createName.localId)
            });
        }

        return completeCreateNamesCreationPayloads;
    }
);

export const selectMagmaIncrementRequestStatusWithCreateNameGroupLocalId = (state: State, createNameGroupLocalId: string): Status | undefined => {
    const requestState = state.names.createNameGroups.magmaIncrementCounterRequestsByLocalId[createNameGroupLocalId];
    if (!requestState) {
        return;
    }
    return requestState.status;
};