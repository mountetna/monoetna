import { CreateName, CreateNameGroup, CreateNameTokenValue } from '../../models';
import { RulesStateSliceForCompleteCreateNames } from '../../selectors/global';
import { NamesState } from './reducer';
import { listToIdObject, listToIdGroupObject, defaultDict } from '../../utils/object';
import { Status } from '../../utils/models';
import { addCreateNameTokenValues, deleteCreateNameTokenValues } from './create-name-token-values';
import { addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds, removeCompleteCreateNamesAndParentsForCreateNameGroupLocalIds } from './complete-create-names';
import { SearchReplaceCriteria, createSearchReplaceCriteriaFromGroups } from '../../utils/names';
import { difference, intersection } from '../../utils/set';



export function addNamesWithGroupsAndTokensValues(
    createNames: CreateName[],
    createNameGroups: CreateNameGroup[],
    createNameTokenValues: CreateNameTokenValue[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {

    const newByRuleName = {...state.createNames.byRuleName};
    createNames.forEach(cn => {
        if (!(cn.ruleName in newByRuleName)) {
            newByRuleName[cn.ruleName] = [];
        }

        newByRuleName[cn.ruleName].push(cn.localId);
    });

    let newState = {
        ...state,
        createNames: {
            byLocalId: {
                ...state.createNames.byLocalId,
                ...listToIdObject(createNames, 'localId'),
            },
            byCreateNameGroupLocalId: {
                ...state.createNames.byCreateNameGroupLocalId,
                ...listToIdGroupObject(createNames, 'createNameGroupLocalId', 'localId'),
            },
            byRuleName: newByRuleName,
        },
        createNameGroups: {
            ...state.createNameGroups,
            byLocalId: {
                ...state.createNameGroups.byLocalId,
                ...listToIdObject(createNameGroups, 'localId'),
            }
        },
    };

    newState.createNameTokenValues = addCreateNameTokenValues(
        createNameTokenValues,
        false,
        rulesStateSliceForCompleteCreateNames,
        newState,
    ).createNameTokenValues;

    const cngLocalIds = createNameGroups.map(cng => cng.localId);

    if (updateCompletionStatus) {
        newState = addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
            cngLocalIds,
            rulesStateSliceForCompleteCreateNames,
            newState,
        );
    }

    return newState;
}

export function deleteGroupsWithNames(
    createNameGroupIds: string[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newGroupsById = { ...state.createNameGroups.byLocalId };
    const newNamesById = { ...state.createNames.byLocalId };
    const newNamesByGroupId = { ...state.createNames.byCreateNameGroupLocalId };
    const newNamesByRuleName = {...state.createNames.byRuleName};

    const cntvLocalIdsToDelete: string[] = [];

    createNameGroupIds.forEach((cngId) => {
        delete newGroupsById[cngId];

        state.createNames.byCreateNameGroupLocalId[cngId].forEach((cnId) => {
            const cn = state.createNames.byLocalId[cnId];
            
            delete newNamesById[cnId];

            // keep track of CreateNameTokenValue.localIds to cleanup
            cntvLocalIdsToDelete.push(...(state.createNameTokenValues.byCreateNameLocalId[cnId] || []));

            newNamesByRuleName[cn.ruleName] = newNamesByRuleName[cn.ruleName].filter(cnLocalId => cnLocalId != cnId);
        });
        delete newNamesByGroupId[cngId];
    });

    let newState = { ...state };

    if (updateCompletionStatus) {
        newState = removeCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(createNameGroupIds, newState);
    }

    return {
        ...newState,
        createNames: {
            byLocalId: newNamesById,
            byCreateNameGroupLocalId: newNamesByGroupId,
            byRuleName: newNamesByRuleName,
        },
        createNameTokenValues: deleteCreateNameTokenValues(
            cntvLocalIdsToDelete,
            false,
            rulesStateSliceForCompleteCreateNames,
            newState,
        ).createNameTokenValues,
        createNameGroups: {
            ...newState.createNameGroups,
            byLocalId: newGroupsById,
            searchLocalIds: removeGroupsFromSearch(createNameGroupIds, newState).createNameGroups.searchLocalIds,
            filterLocalIds: removeGroupsFromFilter(createNameGroupIds, newState).createNameGroups.filterLocalIds,
            replaceLocalIds: removeGroupsFromReplace(createNameGroupIds, newState).createNameGroups.replaceLocalIds,
            selectionLocalIds: removeGroupsFromSelection(createNameGroupIds, newState).createNameGroups.selectionLocalIds,
            magmaCheckDuplicateNameRequestsByLocalId: removeGroupsFromMagmaCheckDuplicateNameRequests(createNameGroupIds, newState).createNameGroups.magmaCheckDuplicateNameRequestsByLocalId,
            magmaIncrementCounterRequestsByLocalId: removeGroupsFromMagmaIncrementCounterRequests(createNameGroupIds, state).createNameGroups.magmaIncrementCounterRequestsByLocalId,
        },
    };
}

export function deleteSelectedGroupsWithNames(
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    return deleteGroupsWithNames(
        [...state.createNameGroups.selectionLocalIds],
        updateCompletionStatus,
        rulesStateSliceForCompleteCreateNames,
        state
    );
}

export function addGroupsToSelection(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            selectionLocalIds: new Set([...state.createNameGroups.selectionLocalIds, ...createNameGroupIds]),
        }
    };
}

export function removeGroupsFromSelection(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            selectionLocalIds: difference(state.createNameGroups.selectionLocalIds, new Set(createNameGroupIds)),
        }
    };
}

export function deselectAllGroups(state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            selectionLocalIds: new Set(),
        }
    };
}

export function setRuleCounterValueForCreateNames(
    ruleCounterValuesByCreateNameLocalId: Record<string, number | undefined>,
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newCreateNamesByLocalId = { ...state.createNames.byLocalId };

    for (const [cnLocalId, ruleCounterValue] of Object.entries(ruleCounterValuesByCreateNameLocalId)) {

        newCreateNamesByLocalId[cnLocalId] = {
            ...newCreateNamesByLocalId[cnLocalId],
            ruleCounterValue
        };
    }

    let newState = {
        ...state,
        createNames: {
            ...state.createNames,
            byLocalId: newCreateNamesByLocalId,
        },
    };

    const cngLocalIds = new Set<string>();

    for (const cnLocalId of Object.keys(ruleCounterValuesByCreateNameLocalId)) {
        cngLocalIds.add(
            state.createNames.byLocalId[cnLocalId].createNameGroupLocalId
        );
    }

    if (updateCompletionStatus === true) {
        newState = addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds(
            [...cngLocalIds],
            rulesStateSliceForCompleteCreateNames,
            newState,
        );
    }

    return newState;
}

// TODO add tests
// TODO add procedure breakdown
export function getMatchedGroupIdsFromSearchCriteria(
    searchCriteria: SearchReplaceCriteria,
    state: NamesState,
    respectFilter: boolean = true
): string[] {

    const searchCriteriaCreateNameTokenValues = Object.values(searchCriteria.byRuleName)
        .map(searchCriteriaByRule => {
            return searchCriteriaByRule.createNameTokenValueLocalIds
                .map(cntvLocalId => state.createNameTokenValues.byLocalId[cntvLocalId]);
        })
        .reduce((cntvList, cntvsToAdd) => [...cntvList, ...cntvsToAdd], []);

    const groupIdsToMatchingValuesCount = defaultDict<string, number>(_ => 0);
    const createNamesToCheckRuleCounter = new Set<CreateName>();

    searchCriteriaCreateNameTokenValues.forEach(searchCriteriaCntv => {
        const cntvIdsToCheck = state.createNameTokenValues.byTokenValueLocalId[searchCriteriaCntv.tokenValueLocalId];

        cntvIdsToCheck.forEach(cntvId => {
            const cntvToCheck = state.createNameTokenValues.byLocalId[cntvId];
            const createName = state.createNames.byLocalId[cntvToCheck.createNameLocalId];
            const ruleSearchCriteria = searchCriteria.byRuleName[createName.ruleName];

            if (
                !ruleSearchCriteria
                || (
                    respectFilter
                    && state.createNameGroups.filterEnabled
                    && !state.createNameGroups.filterLocalIds.has(createName.createNameGroupLocalId)
                )
                || state.createNameGroups.searchLocalIds.has(createName.createNameGroupLocalId)
                || state.createNameGroups.replaceLocalIds.has(createName.createNameGroupLocalId)
            ) {
                return;
            }

            // match if same token value and position in rule
            if (cntvToCheck.ruleTokenLocalId == searchCriteriaCntv.ruleTokenLocalId) {
                ++groupIdsToMatchingValuesCount[createName.createNameGroupLocalId];
            }

            // createNamesToCheckRuleCounter.add(createName);
        });
    });

    // match if same ruleCounterValue
    // (relies on the fact that a Rule will show up only once per CreateNameGroup)
    Object.entries(searchCriteria.byRuleName).forEach(([ruleName, ruleSearchCriteria]) => {
        state.createNames.byRuleName[ruleName].forEach(cnLocalId => {
            const cn = state.createNames.byLocalId[cnLocalId];

            if (
                ruleSearchCriteria.ruleCounterValue != undefined
                && cn.ruleCounterValue == ruleSearchCriteria.ruleCounterValue
            ) {
                ++groupIdsToMatchingValuesCount[cn.createNameGroupLocalId];
            }
        });
    });

    const targetMatchCount = Object.values(searchCriteria.byRuleName)
        .reduce((totalCount, ruleSearchCriteria): number => {

            const ruleCount = ruleSearchCriteria.createNameTokenValueLocalIds.length
                + (ruleSearchCriteria.ruleCounterValue != undefined ? 1 : 0);

            return totalCount + ruleCount;
        }, 0);



    const matchingGroupIds: string[] = [];

    Object.entries(groupIdsToMatchingValuesCount).forEach(([groupId, matchCount]) => {
        if (matchCount == targetMatchCount) {
            matchingGroupIds.push(groupId);
        }
    });

    return matchingGroupIds;
}

export function setGroupsSelectionFromSearchCriteria(state: NamesState): NamesState {
    const searchCriteriaList = createSearchReplaceCriteriaFromGroups(state, state.createNameGroups.searchLocalIds);
    state = deselectAllGroups(state);

    searchCriteriaList.forEach(searchCriteria => {
        state = addGroupsToSelection(getMatchedGroupIdsFromSearchCriteria(searchCriteria, state), state);
    });

    return state;
}

export function setGroupsFilterFromSearchCriteria(state: NamesState): NamesState {
    const searchCriteriaList = createSearchReplaceCriteriaFromGroups(state, state.createNameGroups.searchLocalIds);
    state = disableGroupFilter(state);

    searchCriteriaList.forEach(searchCriteria => {
        state = addGroupsToFilter(getMatchedGroupIdsFromSearchCriteria(searchCriteria, state), state);
    });

    return state;
}

export function addGroupsToSearch(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            searchLocalIds: new Set([...state.createNameGroups.searchLocalIds, ...createNameGroupIds]),
        }
    };
}

export function removeGroupsFromSearch(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            searchLocalIds: difference(state.createNameGroups.searchLocalIds, new Set(createNameGroupIds)),
        }
    };
}

export function addGroupsToReplace(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            replaceLocalIds: new Set([...state.createNameGroups.replaceLocalIds, ...createNameGroupIds]),
        }
    };
}

export function removeGroupsFromReplace(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            replaceLocalIds: difference(state.createNameGroups.replaceLocalIds, new Set(createNameGroupIds)),
        }
    };
}

export function addGroupsToFilter(createNameGroupIds: string[], state: NamesState): NamesState {
    // also need to remove groups from selection if not in filter
    const newFilterIds = new Set([...state.createNameGroups.filterLocalIds, ...createNameGroupIds]);
    const newSelectionIds = intersection(state.createNameGroups.selectionLocalIds, newFilterIds);
    const newSelectionState = addGroupsToSelection([...newSelectionIds], deselectAllGroups(state));

    return {
        ...newSelectionState,
        createNameGroups: {
            ...newSelectionState.createNameGroups,
            filterLocalIds: newFilterIds,
            filterEnabled: true,
        }
    };
}

export function removeGroupsFromFilter(createNameGroupIds: string[], state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            filterLocalIds: difference(state.createNameGroups.filterLocalIds, new Set(createNameGroupIds)),
        }
    };
}

export function disableGroupFilter(state: NamesState): NamesState {
    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            filterLocalIds: new Set(),
            filterEnabled: false,
        }
    };
}

export function setMagmaCheckDuplicateNameRequestForCreateNameGroup(createNameGroupLocalId: string, status: Status, state: NamesState, hasDuplicate?: boolean): NamesState {
    const newComposeErrorState = {
        ...(state.createNameGroups.magmaCheckDuplicateNameRequestsByLocalId[createNameGroupLocalId] || {}),
        status,
        hasDuplicate,
    };

    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            magmaCheckDuplicateNameRequestsByLocalId: {
                ...state.createNameGroups.magmaCheckDuplicateNameRequestsByLocalId,
                [createNameGroupLocalId]: newComposeErrorState,
            },
        },
    };
}

export function removeGroupsFromMagmaCheckDuplicateNameRequests(createNameGroupLocalIds: string[], state: NamesState): NamesState {
    const newmagmaCheckDuplicateNameRequestsByLocalId = { ...state.createNameGroups.magmaCheckDuplicateNameRequestsByLocalId };

    for (const cngLocalId of createNameGroupLocalIds) {
        delete newmagmaCheckDuplicateNameRequestsByLocalId[cngLocalId];
    }

    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            magmaCheckDuplicateNameRequestsByLocalId: newmagmaCheckDuplicateNameRequestsByLocalId,
        },
    };
}

export function setMagmaIncrementCounterRequestForCreateNameGroup(createNameGroupLocalId: string, status: Status, state: NamesState): NamesState {
    const newRequestState = {
        ...(state.createNameGroups.magmaIncrementCounterRequestsByLocalId[createNameGroupLocalId] || {}),
        status,
    };

    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            magmaIncrementCounterRequestsByLocalId: {
                ...state.createNameGroups.magmaIncrementCounterRequestsByLocalId,
                [createNameGroupLocalId]: newRequestState,
            },
        },
    };
}

export function removeGroupsFromMagmaIncrementCounterRequests(createNameGroupLocalIds: string[], state: NamesState): NamesState {
    const newmagmaIncrementCounterRequestsByLocalId = { ...state.createNameGroups.magmaIncrementCounterRequestsByLocalId };

    for (const cngLocalId of createNameGroupLocalIds) {
        delete newmagmaIncrementCounterRequestsByLocalId[cngLocalId];
    }

    return {
        ...state,
        createNameGroups: {
            ...state.createNameGroups,
            magmaIncrementCounterRequestsByLocalId: newmagmaIncrementCounterRequestsByLocalId,
        },
    };
}