import React, { useEffect, useState } from 'react';
import { useSelector, batch } from 'react-redux';
import { makeStyles, Theme } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import ButtonBase from '@material-ui/core/ButtonBase';
import Checkbox from '@material-ui/core/Checkbox';
import Tooltip from '@material-ui/core/Tooltip';
import FileCopyOutlinedIcon from '@material-ui/icons/FileCopyOutlined';
import DeleteOutlineOutlinedIcon from '@material-ui/icons/DeleteOutlineOutlined';
import ErrorOutlineIcon from '@material-ui/icons/ErrorOutline';
import _ from 'lodash';

import { CreateName, CreateNameGroup, CreateNameTokenValue, Rule, RuleToken, TokenValue, UNSET_TOKEN_VALUE, UNSET_VALUE } from '../../../models';
import { selectRulesByName, selectTokenValuesByLocalId, selectTokens, selectRuleTokenLocalIdsWithRuleName, selectRuleTokensByLocalId, selectRuleNamesHierarchicalListByPrimaryRuleName } from '../../../selectors/rules';
import { selectCreateNamesByLocalId, selectCreateNameWithLocalId, selectCreateNameLocalIdsWithGroupId, selectCreateNameTokenValueLocalIdsWithCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectSelectedCreateNameGroupIds, selectSortedCompleteCreateNamesWithCreateNameGroupLocalId, selectCreateNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId, selectRenderedCompleteCreateNamesByLocalId, selectSearchVisible, selectSearchReplaceCriteriaFromSearchGroups, selectReplaceVisible, selectSearchReplaceCriteriaFromReplaceGroups, selectReplaceCreateNameGroupIds } from '../../../selectors/names';
import { addOrReplaceCreateNameTokenValues, setCreateNameRuleCounterValues, duplicateCreateNameGroups, deleteGroupsWithNames, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection, deleteCreateNameTokenValue, setCreateNameGroupComposeError } from '../../../actions/names';
import { selectPathParts } from '../../../selectors/location';
import { TokenSelect } from './select';
import RuleCounterField from './rule-counter-input';
import { createLocalId } from '../../../utils/models';
import { useDispatch } from '../../../utils/redux';
import { selectGlobalState } from '../../../selectors/global';
import { SearchReplaceCriteria, fetchWhetherNameExistsInMagma } from '../../../utils/names';


const useEditorStyles = makeStyles((theme) => ({
    elementsEditorContainer: {
        display: 'inline-flex',
        alignItems: 'center',
        flexWrap: 'nowrap',
        // account for absolute positioning of <RuleCounterField>
        marginRight: '0.4em',
    },
}));


const CreateNameElementsEditor = ({ createName, rule, includeUnsetAsValue, includeRuleCounterIncrementer, parentCompleteCreateNameLocalId, error, highlightMatches }: {
    createName: CreateName,
    rule: Rule,
    includeUnsetAsValue: boolean,
    includeRuleCounterIncrementer: boolean,
    parentCompleteCreateNameLocalId: string | undefined,
    error?: boolean,
    highlightMatches: boolean,
}) => {
    const classes = useEditorStyles();
    const dispatch = useDispatch();

    const globalState = useSelector(selectGlobalState);
    const projectName = useSelector(selectPathParts)[0];
    const tokens = useSelector(selectTokens);

    const ruleTokensByLocalId = useSelector(selectRuleTokensByLocalId);
    let sortedRuleTokens = useSelector(selectRuleTokenLocalIdsWithRuleName(rule.name))
        .map(rtLocalId => ruleTokensByLocalId[rtLocalId]);

    // copying array bc sort modifies in place
    sortedRuleTokens = [...sortedRuleTokens].sort(rt => rt.ord);

    const createNameTokenValuesByLocalId = useSelector(selectCreateNameTokenValuesByLocalId);
    const createNameTokenValues = (useSelector(selectCreateNameTokenValueLocalIdsWithCreateNameLocalId(createName.localId)) || [])
        .map(cntvLocalId => createNameTokenValuesByLocalId[cntvLocalId]);

    const tokenValuesByLocalId = useSelector(selectTokenValuesByLocalId);

    const sortedTokenValues = sortedRuleTokens.map(ruleToken => {
        const createNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        );
        return createNameTokenValue ?
            tokenValuesByLocalId[createNameTokenValue.tokenValueLocalId] : undefined;
    });

    const searchVisible = useSelector(selectSearchVisible);
    const searchCriteria = useSelector(selectSearchReplaceCriteriaFromSearchGroups);
    const replaceVisible = useSelector(selectReplaceVisible);
    const replaceCriteria = useSelector(selectSearchReplaceCriteriaFromReplaceGroups);
    const selectedCreateNameGroupsIds = useSelector(selectSelectedCreateNameGroupIds);
    const replaceCreateNameGroupIds = useSelector(selectReplaceCreateNameGroupIds);

    const renderedTokens: string | undefined = _.every(sortedTokenValues, tv => tv != undefined)
        // @ts-ignore
        ? sortedTokenValues.map(tv => tv.name).reduce((prev, curr) => prev + curr, '')
        : undefined;


    const setTokenValue = (ruleToken: RuleToken, value?: TokenValue) => {
        const oldCreateNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        );

        if (!value) {
            if (oldCreateNameTokenValue) {
                dispatch(deleteCreateNameTokenValue(oldCreateNameTokenValue, globalState));
            }
            return;
        }

        const createNameTokenValue: CreateNameTokenValue = {
            localId: createLocalId(),
            tokenValueLocalId: value.localId,
            createNameLocalId: createName.localId,
            ruleTokenLocalId: ruleToken.localId,
        };

        if (value == UNSET_TOKEN_VALUE) {
            createNameTokenValue.tokenValueLocalId = UNSET_VALUE;
        }

        dispatch(
            addOrReplaceCreateNameTokenValues(
                [createNameTokenValue],
                oldCreateNameTokenValue ? [oldCreateNameTokenValue.localId] : [],
                globalState,
            )
        );
    };

    const setRuleCounterValue = (value?: number) => {
        dispatch(
            setCreateNameRuleCounterValues(
                { [createName.localId]: value },
                globalState,
            )
        );
    };

    const shouldHighlightTokenValue = (ruleToken: RuleToken, tokenValue?: TokenValue): boolean => {
        if (searchVisible) {
            for (const criterium of searchCriteria) {
                const ruleCriterium = criterium.byRuleName[ruleToken.ruleName];

                if (ruleCriterium == undefined) {
                    continue;
                }

                for (const cntvLocalId of ruleCriterium.createNameTokenValueLocalIds) {
                    const cntv = createNameTokenValuesByLocalId[cntvLocalId];

                    if (cntv.tokenValueLocalId == tokenValue?.localId) {
                        return true;
                    }
                }
            }
        }


        if (
            replaceVisible
            && (
                selectedCreateNameGroupsIds.has(createName.createNameGroupLocalId)
                || replaceCreateNameGroupIds.has(createName.createNameGroupLocalId)
            )
        ) {
            for (const criterium of replaceCriteria) {
                const ruleCriterium = criterium.byRuleName[ruleToken.ruleName];

                if (ruleCriterium == undefined) {
                    continue;
                }

                for (const cntvLocalId of ruleCriterium.createNameTokenValueLocalIds) {
                    const cntv = createNameTokenValuesByLocalId[cntvLocalId];

                    if (cntv.ruleTokenLocalId == ruleToken.localId) {
                        return true;
                    }
                }
            }
        }

        return false;
    };

    const shouldHighlightCounterValue = (ruleName: string, counterValue?: number): boolean => {
        if (searchVisible) {
            for (const criterium of searchCriteria) {
                const ruleCriterium = criterium.byRuleName[ruleName];

                if (
                    ruleCriterium != undefined
                    && ruleCriterium.ruleCounterValue != undefined
                    && ruleCriterium.ruleCounterValue == counterValue
                ) {
                    return true;
                }
            }
        }


        if (
            replaceVisible
            && (
                selectedCreateNameGroupsIds.has(createName.createNameGroupLocalId)
                || replaceCreateNameGroupIds.has(createName.createNameGroupLocalId)
            )
        ) {
            for (const criterium of replaceCriteria) {
                const ruleCriterium = criterium.byRuleName[ruleName];

                if (
                    ruleCriterium != undefined
                    && ruleCriterium.ruleCounterValue != undefined
                ) {
                    return true;
                }
            }
        }

        return false;
    };

    return (
        <FormGroup row className={classes.elementsEditorContainer}>
            {
                sortedRuleTokens.map((ruleToken, idx) =>
                    <React.Fragment key={ruleToken.localId}>
                        <TokenSelect
                            token={tokens[ruleToken.tokenName]}
                            value={sortedTokenValues[idx]}
                            onSetTokenValue={value => setTokenValue(ruleToken, value)}
                            includeUnsetAsValue={includeUnsetAsValue}
                            highlight={shouldHighlightTokenValue(ruleToken, sortedTokenValues[idx])}
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
                    includeRuleCounterIncrementer={includeRuleCounterIncrementer}
                    handleSetCounterValue={setRuleCounterValue}
                    highlight={shouldHighlightCounterValue(rule.name, createName.ruleCounterValue)}
                />
            }
        </FormGroup>
    );
};



interface ComposerStylesProps {
    includeRuleCounterIncrementer: boolean
}



const useComposerStyles = makeStyles<Theme, ComposerStylesProps>((theme) => ({
    container: {
        display: 'inline-flex',
        flexWrap: 'wrap',
        '& > *': {
            // account for absolute positioning of RuleCounter Incrementer
            paddingTop: props => props.includeRuleCounterIncrementer ? '1.5em' : 'none',
        },
    },
    toolsContainer: {
        display: 'inline-flex',
        alignItems: 'center',
        marginRight: '1em',

    },
    checkbox: {
        padding: '0',
        paddingRight: '0.25em'
    },
    editorsContainer: {
        display: 'flex',
        alignItems: 'center',
        position: 'relative',
    },
    infoTooltip: {
        color: 'red',
    },
    infoTooltipContent: {
        fontSize: '14px',
    },
    infoTooltipIconContainer: {
        display: 'flex',
        alignItems: 'center',
        '&:hover': {
            cursor: 'help',
        },
    },
}));


const CreateNameGroupComposer = ({
    createNameGroup, className, includeTools = false, includeRuleCounterIncrementer = true,
    includeUnsetAsValue = false, checkForDuplicates = true, highlightMatches = true,
}: {
    createNameGroup: CreateNameGroup,
    className?: string,
    includeTools?: boolean
    includeRuleCounterIncrementer?: boolean,
    includeUnsetAsValue?: boolean,
    checkForDuplicates?: boolean,
    highlightMatches?: boolean,
}) => {
    const dispatch = useDispatch();
    const classes = useComposerStyles({ includeRuleCounterIncrementer });
    const [duplicateTracker, setDuplicateTracker] = useState({ local: 0, remote: 0 });

    const globalState = useSelector(selectGlobalState);
    const projectName = useSelector(selectPathParts)[0];
    const createNameLocalIds = useSelector(selectCreateNameLocalIdsWithGroupId(createNameGroup.localId));
    const createNamesByLocalId = useSelector(selectCreateNamesByLocalId);
    const primaryCreateName = useSelector(selectCreateNameWithLocalId(createNameGroup.primaryCreateNameLocalId));
    const allRules = useSelector(selectRulesByName);
    const orderedRuleNames = useSelector(selectRuleNamesHierarchicalListByPrimaryRuleName)[primaryCreateName.ruleName];
    const selectedCreateNameGroupsIds = useSelector(selectSelectedCreateNameGroupIds);

    const sortedCompleteCreateNamesWithCreateNameGroupLocalId = useSelector(
        selectSortedCompleteCreateNamesWithCreateNameGroupLocalId(createNameGroup.localId)
    );
    const primaryCompleteCreateName = sortedCompleteCreateNamesWithCreateNameGroupLocalId[sortedCompleteCreateNamesWithCreateNameGroupLocalId.length - 1];
    const createNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId = useSelector(selectCreateNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId);
    const renderedCompleteCreateNamesByLocalId = useSelector(selectRenderedCompleteCreateNamesByLocalId);
    const localInstanceCount = (primaryCompleteCreateName != undefined ?
        createNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId[primaryCompleteCreateName.localId] : []
    ).length;

    // check for local duplicates
    useEffect(() => {
        if (checkForDuplicates) {
            batch(() => {
                setDuplicateTracker(tracker => ({ local: localInstanceCount, remote: tracker.remote }));

                if (localInstanceCount > 1 || duplicateTracker.remote) {
                    dispatch(setCreateNameGroupComposeError(createNameGroup.localId, 'idle', true));
                } else {
                    dispatch(setCreateNameGroupComposeError(createNameGroup.localId, 'idle', false));
                }
            });
        }
    }, [localInstanceCount]);

    // check for remote duplicates when completeness changes
    useEffect(() => {
        async function _checkForDuplicates() {
            dispatch(setCreateNameGroupComposeError(createNameGroup.localId, 'inProgress'));

            if (!primaryCompleteCreateName) {
                batch(() => {
                    setDuplicateTracker({ local: 0, remote: 0 });
                    dispatch(setCreateNameGroupComposeError(createNameGroup.localId, 'idle', false));
                });
                return;
            }

            const renderedName = renderedCompleteCreateNamesByLocalId[primaryCompleteCreateName.localId];

            try {
                const remoteDuplicate = await fetchWhetherNameExistsInMagma(
                    projectName,
                    primaryCreateName.ruleName,
                    renderedName
                );

                batch(() => {
                    setDuplicateTracker({ local: localInstanceCount, remote: remoteDuplicate ? 1 : 0 });

                    if (localInstanceCount > 1 || remoteDuplicate) {
                        dispatch(setCreateNameGroupComposeError(createNameGroup.localId, 'idle', true));
                    } else {
                        dispatch(setCreateNameGroupComposeError(createNameGroup.localId, 'idle', false));
                    }
                });
            } catch (err) {
                console.error(`Error determining whether name "${renderedName}" has remote duplicate: ${err}"`);
                dispatch(setCreateNameGroupComposeError(createNameGroup.localId, 'error'));
            }
        }

        if (checkForDuplicates) {
            _checkForDuplicates();
        }
    }, [primaryCompleteCreateName?.localId]);

    const createErrorMessage = () => {
        const errorMsgs: string[] = [];

        if (duplicateTracker.local > 1) {
            errorMsgs.push('locally');
        }
        if (duplicateTracker.remote) {
            errorMsgs.push('in Magma');
        }

        return errorMsgs.length
            ? `Name already exists ${errorMsgs.join(' and ')}`
            : undefined;
    };

    const errorMessage = createErrorMessage();

    const handleClickSelect = (event: React.ChangeEvent<HTMLInputElement>) => {
        if (event.target.checked) {
            dispatch(addCreateNameGroupsToSelection([createNameGroup.localId]));
            return;
        }
        dispatch(removeCreateNameGroupsFromSelection([createNameGroup.localId]));
    };

    const handleClickCopy = () => {
        dispatch(duplicateCreateNameGroups([createNameGroup], globalState));
    };

    const handleClickDelete = () => {
        dispatch(deleteGroupsWithNames([createNameGroup.localId], globalState));
    };

    return (
        <div
            className={`${classes.container} ${errorMessage ? 'hasError' : ''} ${className != undefined ? className : ''}`}
        >
            {
                includeTools &&
                <span className={classes.toolsContainer}>
                    <Checkbox
                        checked={selectedCreateNameGroupsIds.has(createNameGroup.localId)}
                        onChange={handleClickSelect}
                        inputProps={{ 'aria-label': 'Select Name' }}
                        className={classes.checkbox}
                    />
                    <Tooltip title='Copy Name' placement='top'>
                        <ButtonBase
                            onClick={handleClickCopy}
                            aria-label={'Copy Name with Values'}
                            disableRipple
                            disableTouchRipple
                        >
                            <FileCopyOutlinedIcon />
                        </ButtonBase>
                    </Tooltip>
                    <Tooltip title='Delete Name' placement='top'>
                        <ButtonBase
                            onClick={handleClickDelete}
                            aria-label={'Delete Name'}
                            disableRipple
                            disableTouchRipple
                        >
                            <DeleteOutlineOutlinedIcon />
                        </ButtonBase>
                    </Tooltip>
                </span>
            }
            <span className={classes.editorsContainer}>
                {
                    orderedRuleNames.map((ruleName, idx) => {
                        const createNameLocalId = createNameLocalIds.find((cnLocalId) => {
                            return createNamesByLocalId[cnLocalId].ruleName == ruleName;
                        });
                        if (!createNameLocalId) {
                            console.error(`Error creating CreateNameElementsEditor. CreateName with ruleName ${ruleName} not found.`);
                            return;
                        }
                        const rule = allRules[ruleName];
                        const parentCompleteCreateName = sortedCompleteCreateNamesWithCreateNameGroupLocalId[idx - 1];

                        return (
                            <React.Fragment key={rule.name}>
                                <CreateNameElementsEditor
                                    createName={createNamesByLocalId[createNameLocalId]}
                                    rule={rule}
                                    parentCompleteCreateNameLocalId={parentCompleteCreateName?.localId}
                                    includeRuleCounterIncrementer={includeRuleCounterIncrementer}
                                    includeUnsetAsValue={includeUnsetAsValue}
                                    error={errorMessage != undefined}
                                    highlightMatches={highlightMatches}
                                />
                            </React.Fragment>
                        );
                    })
                }
                {errorMessage &&
                    <Tooltip
                        className={classes.infoTooltip}
                        title={
                            <span
                                className={classes.infoTooltipContent}
                            >
                                {errorMessage}
                            </span>
                        }
                    >
                        <span className={classes.infoTooltipIconContainer}>
                            <ErrorOutlineIcon />
                        </span>
                    </Tooltip>}
            </span>
        </div>
    );
};


export default CreateNameGroupComposer;