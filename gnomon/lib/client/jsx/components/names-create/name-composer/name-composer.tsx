import React from "react";
import { useSelector } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import ButtonBase from "@material-ui/core/ButtonBase";
import Checkbox from "@material-ui/core/Checkbox";
import FileCopyOutlinedIcon from '@material-ui/icons/FileCopyOutlined';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import * as _ from "lodash"

import { CreateName, CreateNameGroup, CreateNameTokenValue, Rule, RuleParent, RuleToken, Token, TokenValue, UNSET_TOKEN_VALUE, UNSET_VALUE } from "../../../models";
import { selectRulesByName, selectTokenValuesByLocalId, selectTokens, selectRuleParentLocalIdsByRuleName, selectRuleParentsByLocalId, selectRuleTokenLocalIdsWithRuleName, selectRuleTokensByLocalId, selectRuleNamesHierarchicalListByPrimaryRuleName } from "../../../selectors/rules";
import { selectCreateNamesByLocalId, selectCreateNameWithLocalId, selectCreateNameLocalIdsWithGroupId, selectCreateNameTokenValueLocalIdsWithCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectCreateNameTokenValueLocalIdsByCreateNameLocalId, selectCreateNameLocalIdsByGroupId, selectSelectedCreateNameGroupIds, selectCompleteCreateNameParentsByChildLocalId, selectCompleteCreateNamesByCreateNameLocalId, selectCreateNameCompleteCreateNameLocalIdsByCreateNameLocalId } from "../../../selectors/names";
import { addOrReplaceCreateNameTokenValues, setCreateNameRuleCounterValues, duplicateCreateNameGroups, deleteGroupsWithNames, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection, addOrReplaceCompleteCreateNamesAndParentsForCreateNameGroupLocalIds } from "../../../actions/names";
import { selectPathParts } from "../../../selectors/location"
import { TokenSelect } from "./select";
import RuleCounterField from "./rule-counter-input";
import { createLocalId } from "../../../utils/models";
import { useDispatch } from "../../../utils/redux";
import { State } from "../../../store";
import { selectGlobalState } from "../../../selectors/global";


const useStyles = makeStyles((theme) => ({
    createNameElementsEditor: {
        display: "inline-flex"
    }
}));


const CreateNameElementsEditor = ({ createName, rule, includeUnsetAsValue, completedCreateNameParentLocalId }: {
    createName: CreateName,
    rule: Rule,
    includeUnsetAsValue: boolean,
    completedCreateNameParentLocalId: string | undefined,
}) => {
    const classes = useStyles()
    const dispatch = useDispatch()

    const globalState: State = useSelector(selectGlobalState)
    const projectName: string = useSelector(selectPathParts)[0]
    const tokens: Record<string, Token> = useSelector(selectTokens)

    const ruleTokensByLocalId: Record<string, RuleToken> = useSelector(selectRuleTokensByLocalId)
    let sortedRuleTokens: RuleToken[] = useSelector(selectRuleTokenLocalIdsWithRuleName(rule.name))
        .map(rtLocalId => ruleTokensByLocalId[rtLocalId])

    // copying array bc sort modifies in place
    sortedRuleTokens = [...sortedRuleTokens].sort(rt => rt.ord)

    const createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue> = useSelector(selectCreateNameTokenValuesByLocalId)
    const createNameTokenValues: CreateNameTokenValue[] = (useSelector(selectCreateNameTokenValueLocalIdsWithCreateNameLocalId(createName.localId)) || [])
        .map(cntvLocalId => createNameTokenValuesByLocalId[cntvLocalId])

    const tokenValuesByLocalId: Record<string, TokenValue> = useSelector(selectTokenValuesByLocalId)

    const sortedTokenValues: (TokenValue | undefined)[] = sortedRuleTokens.map(ruleToken => {
        const createNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        )
        return createNameTokenValue ?
            tokenValuesByLocalId[createNameTokenValue.tokenValueLocalId] : undefined
    })

    const renderedTokens: string | undefined = _.every(sortedTokenValues, tv => tv != undefined)
        ? sortedTokenValues.map(tv => tv.name).reduce((prev, curr) => prev + curr, "")
        : undefined


    const setTokenValue = (value: TokenValue, ruleToken: RuleToken) => {
        const createNameTokenValue: CreateNameTokenValue = {
            localId: createLocalId(),
            tokenValueLocalId: value.localId,
            createNameLocalId: createName.localId,
            ruleTokenLocalId: ruleToken.localId,
        }

        if (value == UNSET_TOKEN_VALUE) {
            createNameTokenValue.tokenValueLocalId = UNSET_VALUE
        }

        const oldCreateNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        )
        dispatch(
            addOrReplaceCreateNameTokenValues(
                [createNameTokenValue],
                oldCreateNameTokenValue ? [oldCreateNameTokenValue.localId] : [],
                globalState,
            )
        )
    }

    const setRuleCounterValue = (value?: number) => {
        dispatch(
            setCreateNameRuleCounterValues(
                { [createName.localId]: value },
                globalState,
            )
        )
    }

    return (
        // <span className="create-name-elements-editor">
        <FormGroup row className={classes.createNameElementsEditor}>
            {
                sortedRuleTokens.map((ruleToken, idx) =>

                    <React.Fragment key={ruleToken.localId}>
                        <TokenSelect
                            token={tokens[ruleToken.tokenName]}
                            value={sortedTokenValues[idx]}
                            onSetTokenValue={value => setTokenValue(value, ruleToken)}
                            includeUnsetAsValue={includeUnsetAsValue}
                        ></TokenSelect>
                    </React.Fragment>
                )
            }
            {
                rule.hasCounter &&
                <RuleCounterField
                    value={createName.ruleCounterValue}
                    renderedTokensPrefix={renderedTokens}
                    completeCreateNameParentLocalId={completedCreateNameParentLocalId}
                    projectName={projectName}
                    ruleName={rule.name}
                    handleSetCounterValue={setRuleCounterValue}
                />
            }
        </FormGroup>
    )
}


const CreateNameGroupComposer = ({ createNameGroup, includeTools = false, includeUnsetAsValue = false }: {
    createNameGroup: CreateNameGroup,
    includeTools?: boolean
    includeUnsetAsValue?: boolean
}) => {
    const dispatch = useDispatch()

    const globalState: State = useSelector(selectGlobalState)
    const createNameLocalIds: string[] = useSelector(selectCreateNameLocalIdsWithGroupId(createNameGroup.localId))
    const createNamesByLocalId: Record<string, CreateName> = useSelector(selectCreateNamesByLocalId)
    const primaryCreateName: CreateName = useSelector(selectCreateNameWithLocalId(createNameGroup.primaryCreateNameLocalId))
    const allRules: Record<string, Rule> = useSelector(selectRulesByName)
    const orderedRuleNames: string[] = useSelector(selectRuleNamesHierarchicalListByPrimaryRuleName)[primaryCreateName.ruleName]
    const selectedCreateNameGroupsIds: Set<string> = useSelector(selectSelectedCreateNameGroupIds)
    // const createNameCompleteCreateNameLocalIdsByCreateNameLocalId: Record<string, string> = useSelector(selectCreateNameCompleteCreateNameLocalIdsByCreateNameLocalId)

    // const completeCreateNameParentsByChildLocalId: Record<string, string> = useSelector(selectCompleteCreateNameParentsByChildLocalId)

    const handleClickSelect = (event: React.ChangeEvent) => {
        if (event.target.checked) {
            dispatch(addCreateNameGroupsToSelection([createNameGroup.localId]))
            return
        }
        dispatch(removeCreateNameGroupsFromSelection([createNameGroup.localId]))
    }

    const handleClickCopy = () => {
        dispatch(duplicateCreateNameGroups([createNameGroup], globalState))
    }

    const handleClickDelete = () => {
        dispatch(deleteGroupsWithNames([createNameGroup.localId], globalState))
    }

    return (
        <div className="create-name-group-composer">
            {
                includeTools &&
                <span className="create-name-group-composer-tools">
                    <Checkbox
                        checked={selectedCreateNameGroupsIds.has(createNameGroup.localId)}
                        onChange={handleClickSelect}
                        inputProps={{ 'aria-label': 'Select Name' }}
                    />
                    <ButtonBase
                        onClick={handleClickCopy}
                        aria-label={`Copy Name with Values"`}
                        disableRipple
                        disableTouchRipple
                    >
                        <FileCopyOutlinedIcon />
                    </ButtonBase>
                    <ButtonBase
                        onClick={handleClickDelete}
                        aria-label={`Delete Name"`}
                        disableRipple
                        disableTouchRipple
                    >
                        <DeleteOutlineOutlinedIcon />
                    </ButtonBase>
                </span>
            }
            {
                orderedRuleNames.map((ruleName) => {
                    const createNameLocalId = createNameLocalIds.find((cnLocalId) => {
                        return createNamesByLocalId[cnLocalId].ruleName == ruleName
                    })
                    if (!createNameLocalId) {
                        console.error(`Error creating CreateNameElementsEditor. CreateName with ruleName ${ruleName} not found.`)
                        return
                    }
                    const rule = allRules[ruleName]

                    // const completeCreateNameLocalId: string | undefined = createNameCompleteCreateNameLocalIdsByCreateNameLocalId[createNameLocalId]
                    // const completedCreateNameParentLocalId = completeCreateNameLocalId
                    //     ? completeCreateNameParentsByChildLocalId[completeCreateNameLocalId]
                    //     : undefined

                    return (
                        <React.Fragment key={rule.name}>
                            <CreateNameElementsEditor
                                createName={createNamesByLocalId[createNameLocalId]}
                                rule={rule}
                                // completedCreateNameParentLocalId={completedCreateNameParentLocalId}
                                includeUnsetAsValue={includeUnsetAsValue}
                            />
                        </React.Fragment>
                    )
                })
            }
        </div>
    )
}


export default CreateNameGroupComposer;