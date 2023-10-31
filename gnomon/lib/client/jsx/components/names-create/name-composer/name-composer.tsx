import React from "react";
import { useSelector } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import ButtonBase from "@material-ui/core/ButtonBase";
import Checkbox from "@material-ui/core/Checkbox";
import FileCopyOutlinedIcon from '@material-ui/icons/FileCopyOutlined';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import _ from "lodash"

import { CompleteCreateName, CompleteCreateNameParent, CreateName, CreateNameCompleteCreateName, CreateNameGroup, CreateNameTokenValue, Rule, RuleParent, RuleToken, Token, TokenValue, UNSET_TOKEN_VALUE, UNSET_VALUE } from "../../../models";
import { selectRulesByName, selectTokenValuesByLocalId, selectTokens, selectRuleParentLocalIdsByRuleName, selectRuleParentsByLocalId, selectRuleTokenLocalIdsWithRuleName, selectRuleTokensByLocalId, selectRuleNamesHierarchicalListByPrimaryRuleName } from "../../../selectors/rules";
import { selectCreateNamesByLocalId, selectCreateNameWithLocalId, selectCreateNameLocalIdsWithGroupId, selectCreateNameTokenValueLocalIdsWithCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectSelectedCreateNameGroupIds, selectCreateNameCompleteCreateNameLocalIdsByCreateNameLocalId, selectCreateNameCompleteCreateNamesByLocalId, selectCompleteCreateNameParentLocalIdsByChildLocalId, selectCompleteCreateNameParentsByLocalId, selectSortedCompleteCreateNamesWithCreateNameGroupLocalId } from "../../../selectors/names";
import { addOrReplaceCreateNameTokenValues, setCreateNameRuleCounterValues, duplicateCreateNameGroups, deleteGroupsWithNames, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection } from "../../../actions/names";
import { selectPathParts } from "../../../selectors/location"
import { TokenSelect } from "./select";
import RuleCounterField from "./rule-counter-input";
import { createLocalId } from "../../../utils/models";
import { useDispatch } from "../../../utils/redux";
import { State } from "../../../store";
import { selectGlobalState } from "../../../selectors/global";


const useEditorStyles = makeStyles((theme) => ({
    container: {
        display: "inline-flex"
    }
}));


const CreateNameElementsEditor = ({ createName, rule, includeUnsetAsValue, parentCompleteCreateNameLocalId }: {
    createName: CreateName,
    rule: Rule,
    includeUnsetAsValue: boolean,
    parentCompleteCreateNameLocalId: string | undefined,
}) => {
    const classes = useEditorStyles()
    const dispatch = useDispatch()

    const globalState = useSelector(selectGlobalState)
    const projectName = useSelector(selectPathParts)[0]
    const tokens = useSelector(selectTokens)

    const ruleTokensByLocalId = useSelector(selectRuleTokensByLocalId)
    let sortedRuleTokens = useSelector(selectRuleTokenLocalIdsWithRuleName(rule.name))
        .map(rtLocalId => ruleTokensByLocalId[rtLocalId])

    // copying array bc sort modifies in place
    sortedRuleTokens = [...sortedRuleTokens].sort(rt => rt.ord)

    const createNameTokenValuesByLocalId = useSelector(selectCreateNameTokenValuesByLocalId)
    const createNameTokenValues = (useSelector(selectCreateNameTokenValueLocalIdsWithCreateNameLocalId(createName.localId)) || [])
        .map(cntvLocalId => createNameTokenValuesByLocalId[cntvLocalId])

    const tokenValuesByLocalId = useSelector(selectTokenValuesByLocalId)

    const sortedTokenValues = sortedRuleTokens.map(ruleToken => {
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
        <FormGroup row className={classes.container}>
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
                    parentCompleteCreateNameLocalId={parentCompleteCreateNameLocalId}
                    projectName={projectName}
                    ruleName={rule.name}
                    handleSetCounterValue={setRuleCounterValue}
                />
            }
        </FormGroup>
    )
}


const useComposerStyles = makeStyles((theme) => ({
    container: {
        display: "inline-flex",
    }
}));


const CreateNameGroupComposer = ({ createNameGroup, className, includeTools = false, includeUnsetAsValue = false }: {
    createNameGroup: CreateNameGroup,
    className?: string,
    includeTools?: boolean
    includeUnsetAsValue?: boolean
}) => {
    const dispatch = useDispatch()
    const classes = useComposerStyles()

    const globalState = useSelector(selectGlobalState)
    const createNameLocalIds = useSelector(selectCreateNameLocalIdsWithGroupId(createNameGroup.localId))
    const createNamesByLocalId = useSelector(selectCreateNamesByLocalId)
    const primaryCreateName = useSelector(selectCreateNameWithLocalId(createNameGroup.primaryCreateNameLocalId))
    const allRules = useSelector(selectRulesByName)
    const orderedRuleNames = useSelector(selectRuleNamesHierarchicalListByPrimaryRuleName)[primaryCreateName.ruleName]
    const selectedCreateNameGroupsIds = useSelector(selectSelectedCreateNameGroupIds)

    const sortedCompleteCreateNamesWithCreateNameGroupLocalId = useSelector(
        selectSortedCompleteCreateNamesWithCreateNameGroupLocalId(createNameGroup.localId)
    )

    const handleClickSelect = (event: React.ChangeEvent<HTMLInputElement>) => {
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
        <div className={`${classes.container} ${className != undefined ? className : ""}`}>
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
                orderedRuleNames.map((ruleName, idx) => {
                    const createNameLocalId = createNameLocalIds.find((cnLocalId) => {
                        return createNamesByLocalId[cnLocalId].ruleName == ruleName
                    })
                    if (!createNameLocalId) {
                        console.error(`Error creating CreateNameElementsEditor. CreateName with ruleName ${ruleName} not found.`)
                        return
                    }
                    const rule = allRules[ruleName]
                    const parentCompleteCreateName = sortedCompleteCreateNamesWithCreateNameGroupLocalId[idx - 1]

                    return (
                        <React.Fragment key={rule.name}>
                            <CreateNameElementsEditor
                                createName={createNamesByLocalId[createNameLocalId]}
                                rule={rule}
                                parentCompleteCreateNameLocalId={parentCompleteCreateName?.localId}
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