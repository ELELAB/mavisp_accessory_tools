
from typing import Any, Dict, Iterable, List, Set, Tuple

import networkx as nx

class PathwayGraphFunctions:
    """Utility functions for pathway graph construction and analysis."""

    @staticmethod
    def make_pathway_graph(
        reactions_lists: Iterable[Tuple[str, str, str]],
        pathway_id: str
    ) -> Tuple[nx.DiGraph, List[str], List[str], List[Dict[str, Any]], List[Dict[str, Any]]]:

        G = nx.DiGraph()
        starting_nodes: List[str] = []
        ending_nodes: List[str] = []
        edge_rows: List[Dict[str, Any]] = []

        for current, next_reaction, previous in reactions_lists:

            # Skip invalid nodes to avoid adding empty nodes to the graph.
            if not current:
                print(f"[WARNING] Skipping invalid current: {current}, next={next_reaction}, previous={previous}")
                continue

            if next_reaction:
                G.add_edge(current, next_reaction)
                edge_rows.append({
                    "pathway_id": pathway_id,
                    "source": current,
                    "target": next_reaction
                })

            if previous:
                G.add_edge(previous, current)
                edge_rows.append({
                    "pathway_id": pathway_id,
                    "source": previous,
                    "target": current
                })

            if not previous:
                starting_nodes.append(current)

            if not next_reaction:
                ending_nodes.append(current)

        # Remove duplicated edges while preserving row structure.
        edge_rows = list({(e["pathway_id"], e["source"], e["target"]): e for e in edge_rows}.values())

        node_rows = []
        for node in G.nodes:
            node_rows.append({
                "pathway_id": pathway_id,
                "node": node,
                "is_start": node in starting_nodes,
                "is_end": node in ending_nodes
            })

        return G, starting_nodes, ending_nodes, edge_rows, node_rows

    @staticmethod
    def remove_duplicates_order(lists: Iterable[List[Any]]) -> List[List[Any]]:
        """Remove duplicate lists while preserving order."""
        seen: Set[Tuple[Any, ...]] = set()
        result: List[List[Any]] = []
        for sub in lists:
            t = tuple(sub)
            if t not in seen:
                seen.add(t)
                result.append(list(sub))
        return result

    @staticmethod
    def is_subsequence(sub: List[Any], main: List[Any]) -> bool:
        """Return True if ``sub`` is a subsequence of ``main``."""
        it = iter(main)
        return all(x in it for x in sub)

    @staticmethod
    def find_non_subsequences(lists: Iterable[List[Any]]) -> List[List[Any]]:
        """Return those lists that are not subsequences of any other list."""
        lists = list(lists)
        non_subsequences: List[List[Any]] = []
        for i, sub in enumerate(lists):
            is_sub = False
            for j, main in enumerate(lists):
                if i != j and PathwayGraphFunctions.is_subsequence(sub, main):
                    is_sub = True
                    break
            if not is_sub:
                non_subsequences.append(sub)
        return non_subsequences

