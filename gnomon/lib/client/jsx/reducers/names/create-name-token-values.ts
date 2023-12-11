import { CreateNameTokenValue } from '../../models';
import { RulesStateSliceForCompleteCreateNames } from '../../selectors/global';
import { addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds } from './complete-create-names';
import { NamesState } from './reducer';



export function addCreateNameTokenValues(
    createNameTokenValues: CreateNameTokenValue[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newByCreateNameLocalId = { ...state.createNameTokenValues.byCreateNameLocalId };
    const newByTokenValueLocalId = { ...state.createNameTokenValues.byTokenValueLocalId };
    const newByRuleTokenLocalId = { ...state.createNameTokenValues.byRuleTokenLocalId };
    const newByLocalId = { ...state.createNameTokenValues.byLocalId };

    createNameTokenValues.forEach(cntv => {
        newByLocalId[cntv.localId] = cntv;
        newByCreateNameLocalId[cntv.createNameLocalId] = [
            ...newByCreateNameLocalId[cntv.createNameLocalId] || [],
            cntv.localId
        ];
        newByTokenValueLocalId[cntv.tokenValueLocalId] = [
            ...newByTokenValueLocalId[cntv.tokenValueLocalId] || [],
            cntv.localId
        ];
        newByRuleTokenLocalId[cntv.ruleTokenLocalId] = [
            ...newByRuleTokenLocalId[cntv.ruleTokenLocalId] || [],
            cntv.localId
        ];
    });

    let newState = {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
            byRuleTokenLocalId: newByRuleTokenLocalId,
        }
    };

    const cngLocalIds = new Set<string>();

    for (const cntv of createNameTokenValues) {
        cngLocalIds.add(
            state.createNames.byLocalId[cntv.createNameLocalId].createNameGroupLocalId
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

export function deleteCreateNameTokenValues(
    createNameTokenValueLocalIds: string[],
    updateCompletionStatus: boolean,
    rulesStateSliceForCompleteCreateNames: RulesStateSliceForCompleteCreateNames,
    state: NamesState,
): NamesState {
    const newByCreateNameLocalId = { ...state.createNameTokenValues.byCreateNameLocalId };
    const newByTokenValueLocalId = { ...state.createNameTokenValues.byTokenValueLocalId };
    const newByRuleTokenLocalId = { ...state.createNameTokenValues.byRuleTokenLocalId };
    const newByLocalId = { ...state.createNameTokenValues.byLocalId };

    createNameTokenValueLocalIds.forEach(cntvLocalId => {
        const cntv = state.createNameTokenValues.byLocalId[cntvLocalId];

        delete newByLocalId[cntv.localId];

        if (cntv.createNameLocalId in newByCreateNameLocalId) {
            newByCreateNameLocalId[cntv.createNameLocalId] = newByCreateNameLocalId[cntv.createNameLocalId]
                .filter(cntvLocalId => cntvLocalId != cntv.localId);

            if (newByCreateNameLocalId[cntv.createNameLocalId].length == 0) {
                delete newByCreateNameLocalId[cntv.createNameLocalId];
            }
        }
        if (cntv.tokenValueLocalId in newByTokenValueLocalId) {
            newByTokenValueLocalId[cntv.tokenValueLocalId] = newByTokenValueLocalId[cntv.tokenValueLocalId]
                .filter(cntvLocalId => cntvLocalId != cntv.localId);

            if (newByTokenValueLocalId[cntv.tokenValueLocalId].length == 0) {
                delete newByTokenValueLocalId[cntv.tokenValueLocalId];
            }
        }
        if (cntv.ruleTokenLocalId in newByRuleTokenLocalId) {
            newByRuleTokenLocalId[cntv.ruleTokenLocalId] = newByRuleTokenLocalId[cntv.ruleTokenLocalId]
                .filter(cntvLocalId => cntvLocalId != cntv.localId);

            if (newByRuleTokenLocalId[cntv.ruleTokenLocalId].length == 0) {
                delete newByRuleTokenLocalId[cntv.ruleTokenLocalId];
            }
        }
    });

    let newState = {
        ...state,
        createNameTokenValues: {
            byLocalId: newByLocalId,
            byCreateNameLocalId: newByCreateNameLocalId,
            byTokenValueLocalId: newByTokenValueLocalId,
            byRuleTokenLocalId: newByRuleTokenLocalId,
        }
    };

    const cngLocalIds = new Set<string>();

    for (const cntvLocalId of createNameTokenValueLocalIds) {
        const cnLocalId = state.createNameTokenValues.byLocalId[cntvLocalId].createNameLocalId;

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