import re
import copy
import contextlib
import io
import time
from collections import defaultdict
from typing import Any, Dict, List, Tuple

import requests

try:
    import reactome2py as rc
    from reactome2py import analysis, content, utils
except ImportError:
    rc = None  # type: ignore

try:
    import pybiopax
except ImportError:
    pybiopax = None  # type: ignore

from reactome_pipeline.uniprot_utils import UniprotFunctions


class ReactomeAnalysisFunctions:
    """Functions for extracting information from BioPAX objects."""

    biopax_cache: Dict[str, Any] = {}
    reactome_call_cache: Dict[Any, Any] = {}
    complex_entities_cache: Dict[str, Any] = {}
    disease_link_cache: Dict[str, Any] = {}
    reaction_target_cache: Dict[Tuple[str, str], bool] = {}
    reaction_information_cache: Dict[Tuple[str, Tuple[str, ...]], List[Dict[str, Any]]] = {}

    @staticmethod
    def safe_model_from_reactome_cached(st_id: str, retries: int = 5, sleep: int = 3) -> Any:
        import time
        import requests

        if st_id in ReactomeAnalysisFunctions.biopax_cache:
            return ReactomeAnalysisFunctions.biopax_cache[st_id]

        reactome_id = st_id.split("-")[2]

        for attempt in range(1, retries + 1):
            try:
                model = pybiopax.api.model_from_reactome(reactome_id)
                ReactomeAnalysisFunctions.biopax_cache[st_id] = model
                return model

            except requests.exceptions.RequestException as e:
                print(
                    f"[WARNING] BioPAX download failed for {st_id} "
                    f"({attempt}/{retries}): {e}"
                )

                if attempt < retries:
                    time.sleep(sleep * attempt)

        print(f"[ERROR] Skipping {st_id}: BioPAX model unavailable")
        ReactomeAnalysisFunctions.biopax_cache[st_id] = None
        return None

    @staticmethod
    def safe_reactome_call(func, *args, retries: int = 8, sleep: int = 5, default=None, **kwargs):
        import time
        import requests
        import io
        import contextlib
        import copy

        func_key = f"{getattr(func, '__module__', '')}.{getattr(func, '__name__', str(func))}"

        cache_key = (
            func_key,
            repr(args),
            repr(sorted(kwargs.items()))
        )

        if cache_key in ReactomeAnalysisFunctions.reactome_call_cache:
            return copy.deepcopy(ReactomeAnalysisFunctions.reactome_call_cache[cache_key])

        for attempt in range(1, retries + 1):
            try:
                stdout_buffer = io.StringIO()

                with contextlib.redirect_stdout(stdout_buffer):
                    result = func(*args, **kwargs)

                printed_output = stdout_buffer.getvalue()

                if printed_output:
                    print(printed_output, end="")

                if "Status code returned a value of 404" in printed_output or "Not Found" in printed_output:
                    print(f"[INFO] Reactome returned 404 for {func_key}, args={args}, kwargs={kwargs}. Skipping.")

                    ReactomeAnalysisFunctions.reactome_call_cache[cache_key] = copy.deepcopy(default)
                    return default

                if "Status code returned a value of 429" in printed_output or "Too Many Requests" in printed_output:
                    wait_time = max(60, sleep * attempt * 10)

                    print(
                        f"[WARNING] Reactome rate limit reached for {func_key}, "
                        f"args={args}, kwargs={kwargs} ({attempt}/{retries}). "
                        f"Waiting {wait_time} seconds before retrying."
                    )

                    if attempt < retries:
                        time.sleep(wait_time)
                        continue

                    print(
                        f"[WARNING] Reactome rate limit persisted after {retries} attempts. "
                        "Returning default."
                    )
                    return default


                if (
                    "Status code returned a value of 500" in printed_output
                    or "Internal Server Error" in printed_output
                    or "Status code returned a value of 521" in printed_output
                ):
                    print(
                        f"[WARNING] Reactome server error for {func_key}, args={args}, kwargs={kwargs} "
                        f"({attempt}/{retries}). Retrying."
                    )

                    if attempt < retries:
                        time.sleep(sleep * attempt)
                        continue

                    print(f"[WARNING] Reactome server error persisted after {retries} attempts. Returning default.")
                    return default

                if result is None:
                    print(
                        f"[INFO] Empty response for {func_key}, args={args}, kwargs={kwargs}. "
                        f"No status code detected. Skipping."
                    )
                    return default

                time.sleep(0.5)
                ReactomeAnalysisFunctions.reactome_call_cache[cache_key] = copy.deepcopy(result)
                return result

            except (
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.ReadTimeout,
            ) as e:
                wait_time = max(30, sleep * attempt * 5)
                print(
                    f"[WARNING] Connection problem for {func_key}, "
                    f"args={args}, kwargs={kwargs} ({attempt}/{retries}): {e}. "
                    f"Waiting {wait_time} seconds before retrying."
                )

                if attempt < retries:
                    time.sleep(wait_time)
                    continue

                return default

            except Exception as e:
                error_msg = str(e)

                if "404" in error_msg or "Not Found" in error_msg:
                    print(f"[INFO] Reactome returned 404 for {func_key}, args={args}, kwargs={kwargs}. Skipping.")

                    ReactomeAnalysisFunctions.reactome_call_cache[cache_key] = copy.deepcopy(default)
                    return default
                
                if "429" in error_msg or "Too Many Requests" in error_msg:
                    wait_time = max(60, sleep * attempt * 10)

                    print(
                        f"[WARNING] Reactome rate limit reached for {func_key}, "
                        f"args={args}, kwargs={kwargs} ({attempt}/{retries}). "
                        f"Waiting {wait_time} seconds before retrying."
                    )

                    if attempt < retries:
                        time.sleep(wait_time)
                        continue

                    print(
                        f"[WARNING] Reactome rate limit persisted after {retries} attempts. "
                        "Returning default."
                    )
                    return default

                if "500" in error_msg or "Internal Server Error" in error_msg or "521" in error_msg:
                    print(
                        f"[WARNING] Reactome server error for {func_key}, args={args}, kwargs={kwargs} "
                        f"({attempt}/{retries}). Retrying."
                    )

                    if attempt < retries:
                        time.sleep(sleep * attempt)
                        continue

                    print(f"[WARNING] Reactome server error persisted after {retries} attempts. Returning default.")
                    return default

                print(f"[WARNING] Reactome call failed for {func_key}, args={args}, kwargs={kwargs} ({attempt}/{retries}): {e}")

                if attempt < retries:
                    time.sleep(sleep * attempt)
                    continue

                return default

        return default

    @staticmethod
    def check_attributes(pathway_step_obj: Any) -> None:
        for attr in dir(pathway_step_obj):
            if not attr.startswith("_"):
                try:
                    getattr(pathway_step_obj, attr)
                except Exception:
                    pass
    @staticmethod
    def top_level_pathways(pathways_list: List[Dict[str, Any]]) -> List[Dict[int, Dict[str, str]]]:
        pathways_dic: List[Dict[int, Dict[str, str]]] = []
        for pathway in pathways_list:
            pathway_code = pathway['stId']
            top_level = ReactomeAnalysisFunctions.safe_reactome_call(
                rc.content.event_ancestors,
                id=pathway_code,
                default=[]
            )
            dic: Dict[int, Dict[str, str]] = {}
            for level in top_level:
                for index, item in enumerate(level[::-1]):
                    dic[index] = {'name': item['displayName'], 'id': item['stId']}
                    pathways_dic.append(dic)
        unique_dicts: List[Dict[int, Dict[str, str]]] = []
        seen: Set[frozenset] = set()
        for item in pathways_dic:
            frozen_item = frozenset((k, frozenset(v.items())) for k, v in item.items())
            if frozen_item not in seen:
                seen.add(frozen_item)
                unique_dicts.append(item)
        return unique_dicts

    @staticmethod

    def collect_disease_link(stId: str) -> str:
        if stId in ReactomeAnalysisFunctions.disease_link_cache:
            return ReactomeAnalysisFunctions.disease_link_cache[stId]

        links: List[str] = []

        disease = ReactomeAnalysisFunctions.safe_reactome_call(
            rc.content.query_id,
            id=stId,
            enhanced=True,
            default=None
        )

        if not disease or not disease.get("isInDisease"):
            ReactomeAnalysisFunctions.disease_link_cache[stId] = ""
            return ""

        for x in disease.get("crossReference", []):
            if x.get("databaseName") == "Mondo" and x.get("url"):
                links.append(x.get("url"))

        result = "_".join(links)

        ReactomeAnalysisFunctions.disease_link_cache[stId] = result
        return result

    @staticmethod
    def get_pathway_order(biopax_obj: Any) -> Dict[str, List[Dict[str, List[str]]]]:
        step_processes: Dict[str, List[Dict[str, List[str]]]] = {}
        for pathway_step_obj in biopax_obj.get_objects_by_type(pybiopax.biopax.BiochemicalReaction):
            Step_Process_Of = pathway_step_obj.step_process_of
            single_step_processes: List[Dict[str, List[str]]] = []
            if Step_Process_Of:
                for process in list(Step_Process_Of):
                    single_step_process: Dict[str, List[str]] = {}
                    next_step = process.next_step
                    next_step_of = process.next_step_of
                    if next_step:
                        next_steps: List[str] = []
                        for step in list(next_step):
                            next_steps.extend(i.name for i in step.step_process)
                        single_step_process['next_step'] = next_steps
                    if next_step_of:
                        next_steps_of: List[str] = []
                        for step in list(next_step_of):
                            next_steps_of.extend(i.name for i in step.step_process)
                        single_step_process['previous'] = next_steps_of
                    single_step_processes.append(single_step_process)
            reaction_name = (
                getattr(pathway_step_obj, "display_name", None)
                or " ".join(getattr(pathway_step_obj, "name", []) or [])
            )

            if not reaction_name:
                print("[WARNING] Reaction without display_name/name")
                continue

            step_processes[reaction_name] = single_step_processes
        return step_processes

    @staticmethod
    def biopax_model_contains_uniprot(biopax_model: Any, uniprot_ac: str) -> bool:
        for protein in biopax_model.get_objects_by_type(pybiopax.biopax.Protein):

            for xref in getattr(protein, "xref", []):
                if str(getattr(xref, "id", "")) == uniprot_ac:
                    return True

            entity_reference = getattr(protein, "entity_reference", None)
            if entity_reference is not None:
                for xref in getattr(entity_reference, "xref", []):
                    if str(getattr(xref, "id", "")) == uniprot_ac:
                        return True

        return False

    @staticmethod
    def build_reaction_information(
        reaction: Dict[str, Any],
        reaction_biopax_model: Any,
        pathways_to_not_analyze: List[str]
    ) -> List[Dict[str, Any]]:
        reaction_stID = reaction["stId"]
        reaction_name = reaction["displayName"]

        biochemical_reactions = list(
            reaction_biopax_model.get_objects_by_type(pybiopax.biopax.BiochemicalReaction)
        )

        reaction_components: List[Any] = [
            biochemical_reaction.name
            for biochemical_reaction in biochemical_reactions
        ]

        should_analyze = True

        if len(reaction_components) > 1:
            pathways_component: List[str] = []

            for biochemical_reaction in biochemical_reactions:
                for path_name in list(getattr(biochemical_reaction, "pathway_component_of", [])):
                    pathways_component.append(path_name.display_name)

            if any(
                any(pathway_to_skip in pathway_component for pathway_to_skip in pathways_to_not_analyze)
                for pathway_component in pathways_component
            ):
                should_analyze = False

        if not should_analyze:
            print(
                f"[INFO] Skipping {reaction_stID}: contains pathway components marked as not to analyze."
            )
            return []

        biochemical_information = ReactomeAnalysisFunctions.get_biochemical(
            reaction_biopax_model
        )
        complexes_information = ReactomeAnalysisFunctions.get_complex(
            reaction_biopax_model
        )
        control_information = ReactomeAnalysisFunctions.get_control_information(
            reaction_biopax_model
        )
        pathway_information = ReactomeAnalysisFunctions.get_pathway(
            reaction_biopax_model
        )
        proteins_information = ReactomeAnalysisFunctions.get_protein(
            reaction_biopax_model
        )
        disease = ReactomeAnalysisFunctions.collect_disease_link(
            reaction_stID
        )

        if len(reaction_components) > 1:
            display_names = reaction_components
        else:
            display_names = [reaction_name]

        reaction_information_list: List[Dict[str, Any]] = []

        for reaction_name_component in display_names:
            if isinstance(reaction_name_component, list):
                display_name = " ".join(map(str, reaction_name_component))
            else:
                display_name = str(reaction_name_component)

            reaction_information: Dict[str, Any] = {
                "stID": reaction_stID,
                "Display_Name": display_name,
                "biochemical": biochemical_information,
                "control_information": control_information,
                #"catalytic": catalysis_information,
                "complex": complexes_information,
                "pathway": pathway_information,
                "proteins": proteins_information,
                "disease": disease
            }

            reaction_information_list.append(reaction_information)

        return reaction_information_list

    @staticmethod
    def get_complex(biopax_obj: Any) -> List[Dict[str, Any]]:
        reaction_complexes: List[Dict[str, Any]] = []
        for pathway_step_obj in biopax_obj.get_objects_by_type(pybiopax.biopax.Complex):
            complexes: Dict[str, Any] = {}
            Name = pathway_step_obj.name
            complex_name = "/".join(Name)
            Component_stoichiometry = pathway_step_obj.component_stoichiometry
            components: List[Dict[str, Any]] = []
            complex_ids: List[str] = []
            for reac_id in pathway_step_obj.xref:
                if reac_id.id.startswith("R-HSA"):
                    complex_ids.append(reac_id.id)
            complex_ids = list(set(complex_ids))
            for stoc in Component_stoichiometry:
                st_id = None
                stoc_coefficient = stoc.stoichiometric_coefficient
                physical_entity = stoc.physical_entity

                for x in getattr(stoc.physical_entity, "xref", []):
                    xref_id = str(getattr(x, "id", ""))

                    if (
                        xref_id.startswith("R-HSA-")
                        or xref_id.startswith("R-ALL-")
                        or xref_id.startswith("R-COV-")
                        or xref_id.startswith("R-HIV-")
                        or xref_id.startswith("R-FLU-")
                        or xref_id.startswith("R-ECO-")
                        or xref_id.startswith("R-MTU-")
                    ):
                        st_id = xref_id
                        break

                components.append({
                    "physical_entity": physical_entity,
                    "stoc_coefficient": stoc_coefficient,
                    "stId": st_id
                })
            complexes["complex_name"] = complex_name
            complexes["complexes_ids"] = complex_ids
            complexes["complex_components"] = components
            Participant = pathway_step_obj.participant_of
            complexes['Participant'] = [i.display_name for i in Participant]
            Evidence = pathway_step_obj.evidence
            complexes['evidence'] = Evidence
            reaction_complexes.append(complexes)
        return reaction_complexes

    @staticmethod
    def get_protein(biopax_obj: Any) -> List[Dict[str, Any]]:
        proteins: List[Dict[str, Any]] = []
        for pathway_step_obj in biopax_obj.get_objects_by_type(pybiopax.biopax.Protein):
            protein_features: Dict[str, Any] = {}
            cellular_location = pathway_step_obj.cellular_location
            entity_reference = pathway_step_obj.entity_reference
            component_of_info = []

            for comp in pathway_step_obj.component_of:
                comp_id = None

                for ref in getattr(comp, "xref", []):
                    ref_id = getattr(ref, "id", None)
                    if ref_id and ref_id.startswith("R-HSA"):
                        comp_id = ref_id
                        break  # Use the first Reactome identifier found for this component.

                component_of_info.append({
                    "name": getattr(comp, "display_name", str(comp)),
                    "stId": comp_id
                })
            controller_of = pathway_step_obj.controller_of
            display_name = pathway_step_obj.display_name
            evidence = pathway_step_obj.evidence
            feature = pathway_step_obj.feature
            get_plain_names = pathway_step_obj.get_plain_names
            list_types = pathway_step_obj.list_types
            member_physical_entity = pathway_step_obj.member_physical_entity
            member_physical_entity_of = pathway_step_obj.member_physical_entity_of
            names = pathway_step_obj.name
            refs = pathway_step_obj.xref
            stId = None
            for ref in refs:
                if ref.id.startswith("R-HSA"):
                    stId = ref.id
            uniprot_acs: List[str] = []
            if '-' not in display_name:
                if entity_reference:
                    uniprot_acs.append(entity_reference.name[0].split(" ")[0].split(":")[1])
            else:
                splitted_names = display_name.split("-")
                for splitted_name in splitted_names:
                    if splitted_name not in ["p", "S", "T", "Y"]:
                        uniprot_ac = UniprotFunctions.uniprot_gene_to_uniprot_ac(splitted_name)
                        if uniprot_ac and uniprot_ac not in uniprot_acs:
                            uniprot_acs.append(uniprot_ac)
                if len(uniprot_acs) == 0:
                    uniprot_acs = []
                    if entity_reference:
                        uniprot_acs.append(entity_reference.name[0].split(" ")[0].split(":")[1])
            not_feature = pathway_step_obj.not_feature
            participant_of = pathway_step_obj.participant_of
            standard_name = pathway_step_obj.standard_name
            uid = pathway_step_obj.uid
            uid_index = re.search(r'\d+$', uid).group()
            if feature:
                features: List[Dict[str, Any]] = []
                for i in feature:
                    features_dic: Dict[str, Any] = {}
                    features_dic["feature_location"] = i.feature_location
                    features_dic["feature_location_type"] = i.feature_location_type
                    features_dic["feature_of"] = i.feature_of
                    features_dic["uid"] = i.uid
                    features.append(features_dic)
                    if "modification_type" in dir(i):
                        features_dic["modification_type"] = i.modification_type
                protein_features["feature"] = features
            protein_features["cellular_location"] = cellular_location
            protein_features["component_of"] = component_of_info
            protein_features["controller_of"] = controller_of
            protein_features["member_physical_entity"] = member_physical_entity
            protein_features["member_physical_entity_of"] = member_physical_entity_of
            protein_features["display_name"] = display_name
            protein_features["available_names"] = "_".join(names)
            protein_features["stId"] = stId
            protein_features["uniprot_ac"] = "_".join(uniprot_acs)
            proteins.append(protein_features)
        return proteins

    @staticmethod
    def get_control_information(biopax_obj: Any) -> List[Dict[str, Any]]:
        reaction_controllers: List[Dict[str, Any]] = []
        for pathway_step_obj in biopax_obj.get_objects_by_type(pybiopax.biopax.Control):
            reaction_controller: Dict[str, Any] = {}
            control_type = pathway_step_obj.control_type
            controlled = pathway_step_obj.controlled.display_name
            controlled_of = pathway_step_obj.controlled_of
            controller = pathway_step_obj.controller
            controller_stid = []
            for controller_obj in pathway_step_obj.controller:
                for x in controller_obj.xref:
                    if str(x.id).startswith("R-HSA-"):
                        controller_stid.append(x.id)
            display_name = pathway_step_obj.display_name
            get_plain_names = pathway_step_obj.get_plain_names()
            interaction_type = pathway_step_obj.interaction_type
            participant = pathway_step_obj.participant
            participant_of = pathway_step_obj.participant_of
            pathway_component_of = pathway_step_obj.pathway_component_of
            standard_name = pathway_step_obj.standard_name
            evidence = pathway_step_obj.evidence
            step_process_of = pathway_step_obj.step_process_of
            reaction_controller['control_type'] = control_type
            reaction_controller['controlled'] = controlled
            reaction_controller['controlled_of'] = controlled_of
            reaction_controller['controller'] = controller
            reaction_controller['controller_stid'] = controller_stid
            reaction_controller['participant'] = participant
            reaction_controller['pathway_component_of'] = pathway_component_of
            reaction_controller['evidence'] = evidence
            reaction_controller['control_type'] = control_type
            reaction_controller['interaction_type'] = interaction_type
            reaction_controllers.append(reaction_controller)
        return reaction_controllers

    @staticmethod
    def get_catalysis(biopax_obj: Any) -> Dict[str, List[Any]]:
        all_catalysis: Dict[str, List[Any]] = defaultdict(list)

        for pathway_step_obj in biopax_obj.get_objects_by_type(pybiopax.biopax.Catalysis):
            catalysis_direction = pathway_step_obj.catalysis_direction
            cofactor = pathway_step_obj.cofactor
            control_type = pathway_step_obj.control_type
            controlled = pathway_step_obj.controlled.display_name
            controlled_of = pathway_step_obj.controlled_of
            controller = pathway_step_obj.controller
            display_name = pathway_step_obj.display_name
            interaction_type = pathway_step_obj.interaction_type
            participant = pathway_step_obj.participant
            pathway_component_of = pathway_step_obj.pathway_component_of
            standard_name = pathway_step_obj.standard_name
            evidence = pathway_step_obj.evidence

            all_catalysis["display_name"].append(display_name)
            all_catalysis["standard_name"].append(standard_name)
            all_catalysis["catalysis_direction"].append(catalysis_direction)
            all_catalysis["cofactor"].append(cofactor)
            all_catalysis["control_type"].append(control_type)
            all_catalysis["controlled"].append(controlled)
            all_catalysis["controlled_of"].append(controlled_of)
            all_catalysis["controller"].append(controller)
            all_catalysis["participant"].append(participant)
            all_catalysis["pathway_component_of"].append(pathway_component_of)
            all_catalysis["evidence"].append(evidence)
            all_catalysis["interaction_type"].append(interaction_type)

        return dict(all_catalysis)

    @staticmethod
    def get_pathway(biopax_obj: Any) -> None:
        # This function in the original script prints debug information and
        # returns nothing.  It has been retained for completeness but does not
        # contribute to the data collection.
        for pathway_step_obj in biopax_obj.get_objects_by_type(pybiopax.biopax.Pathway):
            dir(pathway_step_obj)
        return None

    @staticmethod
    def get_biochemical(biopax_obj: Any) -> List[Dict[str, Any]]:
        all_reactions_parameters: List[Dict[str, Any]] = []
        for pathway_step_obj in biopax_obj.get_objects_by_type(pybiopax.biopax.BiochemicalReaction):
            reaction_parameters: Dict[str, Any] = {}
            Name = pathway_step_obj.name
            reaction_name = " ".join(Name)
            Display_Name = pathway_step_obj.display_name
            Left = pathway_step_obj.left
            Right = pathway_step_obj.right
            Participant = pathway_step_obj.participant
            Evidence = pathway_step_obj.evidence
            Standard_Name = pathway_step_obj.standard_name
            Conversion_Direction = pathway_step_obj.conversion_direction
            Delta_G = pathway_step_obj.delta_g
            Delta_H = pathway_step_obj.delta_h
            Delta_S = pathway_step_obj.delta_s
            K_eq = pathway_step_obj.k_e_q
            Pathway_Component_Of = pathway_step_obj.pathway_component_of
            pathways_component: List[str] = []
            if Pathway_Component_Of:
                for path_name in list(Pathway_Component_Of):
                    pathways_component.append(path_name.display_name)
            reaction_parameters['Left'] = Left
            reaction_parameters['Right'] = Right
            reaction_parameters['Participant'] = Participant
            reaction_parameters['Pathway_Component_Of'] = pathways_component
            reaction_parameters['Pathways_component'] = pathways_component
            reaction_parameters['Evidence'] = Evidence
            reaction_parameters['Standard_Name'] = Standard_Name
            reaction_parameters['Conversion_Direction'] = Conversion_Direction
            reaction_parameters['K_eq'] = K_eq
            reaction_parameters['Delta_G'] = Delta_G
            reaction_parameters['Delta_H'] = Delta_H
            reaction_parameters['Delta_S'] = Delta_S
            all_reactions_parameters.append(reaction_parameters)
        return all_reactions_parameters

