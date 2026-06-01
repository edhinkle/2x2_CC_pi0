"""
Histogram loader for reading ROOT files and organizing histograms.

Provides functionality to load histograms from ROOT files and organize them
by type (Reco, Truth, TruthMatched, CutFlow) with caching and metadata.
"""

from typing import Dict, List, Optional, Any
from pathlib import Path
import warnings


class HistogramLoader:
    """
    Load and manage histograms from ROOT files.
    
    This class handles:
    - Loading ROOT files using uproot (no ROOT installation needed)
    - Organizing histograms by category (Reco, Truth, Matched, CutFlow)
    - Caching histograms in memory
    - Validating histogram names
    
    Assumes histograms are organized in ROOT directories:
    - root_file["reco"] - Reco histograms
    - root_file["truth"] - Truth histograms
    - root_file["truthMatched"] - TruthMatched histograms
    - root_file["cutflows"] - CutFlow histograms
    
    Attributes:
        file_path (Path): Path to ROOT file
        root_file: Opened ROOT file object
        hist_dict (Dict): Cache of loaded histograms
    """
    
    # Directory names in ROOT file
    RECO_DIR = "reco"
    TRUTH_DIR = "truth"
    MATCHED_DIR = "truthMatched"
    CUTFLOW_DIR = "cutflows"
    
    def __init__(self, file_path: str | Path):
        """
        Initialize loader with a ROOT file.
        
        Args:
            file_path: Path to ROOT file
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ImportError: If uproot is not installed
        """
        try:
            import uproot
        except ImportError:
            raise ImportError(
                "uproot package required. Install with: pip install uproot"
            )
        
        self.file_path = Path(file_path)
        
        if not self.file_path.exists():
            raise FileNotFoundError(f"ROOT file not found: {self.file_path}")
        
        self.root_file = uproot.open(self.file_path)
        self.hist_dict: Dict[str, Any] = {}
    
    def list_histograms(self, directory: Optional[str] = None) -> List[str]:
        """
        List all histogram names in the file or in a specific directory.
        
        Args:
            directory: If specified, list histograms only in that directory
                      (e.g., "reco", "truth", "truthMatched", "cutflows")
        
        Returns:
            List of histogram names
        """
        if directory:
            if directory not in self.root_file:
                raise KeyError(f"Directory '{directory}' not found in ROOT file")
            
            dir_obj = self.root_file[directory]
            # Get all keys in directory, filtering out subdirectories
            return [key for key in dir_obj.keys() if not key.startswith("_")]
        
        return list(self.root_file.keys())
    
    def _get_from_directory(self, directory: str) -> Dict[str, Any]:
        """
        Load all histograms from a directory in the ROOT file.
        
        Args:
            directory: Directory name (e.g., "reco", "truth", "matched", "cutflows")
        
        Returns:
            Dictionary mapping histogram names to histogram objects
            
        Raises:
            KeyError: If directory not found
        """
        if directory not in self.root_file:
            raise KeyError(
                f"Directory '{directory}' not found. "
                f"Available: {list(self.root_file.keys())}"
            )
        
        dir_obj = self.root_file[directory]
        loaded = {}
        
        for name in self.list_histograms(directory):
            try:
                hist_path = f"{directory}/{name}"
                
                # Check cache first
                if hist_path not in self.hist_dict:
                    self.hist_dict[hist_path] = dir_obj[name]
                
                loaded[name] = self.hist_dict[hist_path]
            except Exception as e:
                warnings.warn(f"Failed to load {directory}/{name}: {e}")
        
        return loaded
    
    def get_histogram(self, hist_name: str, directory: Optional[str] = None,
                     use_cache: bool = True) -> Any:
        """
        Load a single histogram by name.
        
        Args:
            hist_name: Name of histogram to load
            directory: Optional directory name (e.g., "reco", "truth", "matched", "cutflows")
            use_cache: Whether to use cached version if available
        
        Returns:
            Histogram object
            
        Raises:
            KeyError: If histogram not found
        """
        if directory:
            hist_path = f"{directory}/{hist_name}"
        else:
            hist_path = hist_name
        
        # Check cache first
        if use_cache and hist_path in self.hist_dict:
            return self.hist_dict[hist_path]
        
        # Load from file
        try:
            if directory:
                if directory not in self.root_file:
                    raise KeyError(f"Directory '{directory}' not found")
                hist = self.root_file[directory][hist_name]
            else:
                hist = self.root_file[hist_name]
            
            self.hist_dict[hist_path] = hist
            return hist
        
        except KeyError:
            if directory:
                available = self.list_histograms(directory)
                raise KeyError(
                    f"Histogram '{hist_name}' not found in '{directory}' directory.\n"
                    f"Available: {available}"
                )
            else:
                available = self.list_histograms()
                raise KeyError(
                    f"Histogram '{hist_name}' not found.\n"
                    f"Available: {available}"
                )
    
    def load_histograms(self, hist_names: List[str], 
                       directory: Optional[str] = None,
                       use_cache: bool = True) -> Dict[str, Any]:
        """
        Load multiple histograms by name.
        
        Args:
            hist_names: List of histogram names
            directory: Optional directory name (e.g., "reco", "truth", "matched", "cutflows")
            use_cache: Whether to use cached versions if available
        
        Returns:
            Dictionary mapping names to histogram objects
        """
        loaded = {}
        for name in hist_names:
            try:
                loaded[name] = self.get_histogram(name, directory=directory,
                                                 use_cache=use_cache)
            except KeyError as e:
                warnings.warn(f"Failed to load {name}: {e}")
        
        return loaded
    
    def get_reco_histograms(self) -> Dict[str, Any]:
        """
        Load all RecoHists histograms from "reco" directory.
        
        Returns:
            Dictionary of RecoHists histograms (histogram name -> histogram object)
        """
        return self._get_from_directory(self.RECO_DIR)
    
    def get_truth_histograms(self) -> Dict[str, Any]:
        """
        Load all TruthHists histograms from "truth" directory.
        
        Returns:
            Dictionary of TruthHists histograms (histogram name -> histogram object)
        """
        return self._get_from_directory(self.TRUTH_DIR)
    
    def get_matched_histograms(self) -> Dict[str, Any]:
        """
        Load all TruthMatchedHists histograms from "matched" directory.
        
        Returns:
            Dictionary of TruthMatchedHists histograms (histogram name -> histogram object)
        """
        return self._get_from_directory(self.MATCHED_DIR)
    
    def get_cutflow_histograms(self) -> Dict[str, Any]:
        """
        Load all CutFlow histograms from "cutflows" directory.
        
        Returns:
            Dictionary of CutFlow histograms (histogram name -> histogram object)
        """
        return self._get_from_directory(self.CUTFLOW_DIR)
    
    def get_all_histograms(self) -> Dict[str, Dict[str, Any]]:
        """
        Load all histograms organized by category.
        
        Returns:
            Dictionary with keys "reco", "truth", "matched", "cutflows", each containing
            a dictionary of histograms
        """
        all_hists = {
            "reco": self.get_reco_histograms(),
            "truth": self.get_truth_histograms(),
            "matched": self.get_matched_histograms(),
            "cutflows": self.get_cutflow_histograms()
        }
        
        return all_hists
    
    def filter_histograms(self, pattern: str,
                         directory: Optional[str] = None) -> Dict[str, Any]:
        """
        Filter histograms by name pattern.
        
        Args:
            pattern: Pattern to match (e.g., "CosL", "Multiplicity")
                    Matching is case-insensitive
            directory: Optional directory to search (e.g., "reco", "truth", "matched", "cutflows")
        
        Returns:
            Dictionary of matching histograms
        """
        pattern_lower = pattern.lower()
        matching = {}
        
        if directory:
            hist_names = self.list_histograms(directory)
            for name in hist_names:
                if pattern_lower in name.lower():
                    matching[name] = self.get_histogram(name, directory=directory)
        else:
            # Search across all directories
            dirs_to_search = [self.RECO_DIR, self.TRUTH_DIR, self.MATCHED_DIR, self.CUTFLOW_DIR]
            for directory_name in dirs_to_search:
                try:
                    hist_names = self.list_histograms(directory_name)
                    for name in hist_names:
                        if pattern_lower in name.lower():
                            key = f"{directory_name}/{name}"
                            matching[key] = self.get_histogram(name, directory=directory_name)
                except KeyError:
                    pass  # Directory doesn't exist, skip
        
        return matching
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about loaded histograms.
        
        Returns:
            Dictionary with statistics
        """
        stats = {
            "file": str(self.file_path),
            "cached_histograms": len(self.hist_dict),
        }
        
        # Count histograms in each directory
        for dir_name in [self.RECO_DIR, self.TRUTH_DIR, self.MATCHED_DIR, self.CUTFLOW_DIR]:
            try:
                count = len(self.list_histograms(dir_name))
                stats[f"{dir_name}_histograms"] = count
            except KeyError:
                pass
        
        return stats
    
    def list_directories(self) -> List[str]:
        """
        List all top-level directories in the ROOT file.
        
        Returns:
            List of directory names
        """
        return [key for key in self.root_file.keys() if not key.startswith("_")]
    
    def close(self) -> None:
        """Close the ROOT file and clear cache."""
        if self.root_file:
            self.root_file.close()
        self.hist_dict.clear()
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
    
    def __repr__(self) -> str:
        """String representation."""
        stats = self.get_statistics()
        return (
            f"HistogramLoader('{self.file_path.name}', "
            f"reco={stats.get('reco_histograms', 0)}, "
            f"truth={stats.get('truth_histograms', 0)}, "
            f"matched={stats.get('matched_histograms', 0)}, "
            f"cutflows={stats.get('cutflows_histograms', 0)})"
        )
